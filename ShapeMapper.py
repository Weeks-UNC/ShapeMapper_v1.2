#!/usr/bin/env python
#  Primary script to run the SHAPE-MaP analysis pipeline (in other words RUN THIS SCRIPT)
#
#  - Requires a config file as an argument (an EXAMPLE.cfg should have come with this script)
#  - Requires the following programs be executable from any location (i.e. in the PATH):
#        python, bowtie2, bowtie2-build
#  - Take a look at the README for other required modules, installation, and execution help.
#  - Public release 1.2
#  - Copyright Steven Busan 2014

##################################################################################
# GPL statement:
# This file is part of Shapemapper.
#
# ShapeMapper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ShapeMapper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ShapeMapper.  If not, see <http://www.gnu.org/licenses/>.
##################################################################################

import sys, os, traceback, time
import subprocess as sp
import parseConfigFile
import pivotCSV

logFilePath = os.path.join(os.getcwd(),"log.txt")
# create master log file for this run
if len(sys.argv) > 2:
    logFilePath = sys.argv[2]
    print "creating log file "+logFilePath
    log = open(logFilePath, "w") # allow the user to specify where the log file should go
                                 # This allows for running multiple pipeline instances in the same dir
                                 # - if doing this, requires the following workarounds (not included
                                 #   in this package):
                                 # 1. Config file pickling creates file access collisions, so start 
                                 #    parallel job slowly one at a time.
                                 # 2. Bowtie2 fails for concurrent access to the same reference index,
                                 #    so create duplicate .fa sequence files with different names for
                                 #    each parallel job.
else:
    print "no log specified, creating default instead"
    log = open(logFilePath,"w")
print "Run details will be written to %s"%logFilePath 

# create directory for temporary stdout, stderr storage
tempDir = os.path.abspath(os.path.join(os.getcwd(),"temp"))
#print "tempDir = "+tempDir
if not os.path.exists(tempDir):
    os.makedirs(tempDir)


def timeStamp():
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return "at %i:%i:%i, %i/%i"%(hour, minute, second, month, day)

log.write("Pipeline started in %s %s\n"%(os.getcwd(),timeStamp()))
log.flush()

def chunk(elements, chunkSize):
    # break iterator into chunks
    while len(elements) >= chunkSize:
        yield elements[:chunkSize]
        elements = elements[chunkSize:]
    if len(elements):
        yield elements

def spawnProcesses(maxProcesses, processNames, argListDict):
# spawn a given number of child processes (spawning more as children complete)
# and return stdouts and stderrs
    processNamesLocal = list(processNames)
    #[os.path.basename(name) for name in list(processNames)]
    processes = {}                                                                                     
    stdOuts = {}
    stdErrs = {}
    stdOutFiles = {}
    stdErrFiles = {}
    while len(processNamesLocal) != 0 or len(processes) != 0:
        time.sleep(1) # make this a less-busy wait         
        if len(processes) < maxProcesses and len(processNamesLocal) > 0:
            procName = processNamesLocal.pop()
            #print "procName = "+procName
            #print "argList = "+str(argListDict[procName])
            stdOutFiles[procName] = open(os.path.join(tempDir,os.path.basename(procName)+".stdout"),"w")
            #print "Created temp file %s for writing."%(os.path.join(tempDir,os.path.basename(procName)+".stdout"))
            stdErrFiles[procName] = open(os.path.join(tempDir,os.path.basename(procName)+".stderr"),"w") 
            # TODO: add timestamp to each stdout as it completes, instead of logging time that all processes complete?
            # or maybe define an output handling callback function and call it after each process completes
            processes[procName] = sp.Popen(argListDict[procName], stdout=stdOutFiles[procName], stderr=stdErrFiles[procName])
        if len(processes.keys()) != 0:
            for procName in processes.keys():
                if processes[procName].poll() is not None:
                    processes.pop(procName)
    for procName in stdOutFiles.keys():
        stdOutFiles[procName].close()
        stdErrFiles[procName].close()
        stdOutFile = open(os.path.join(tempDir,os.path.basename(procName)+".stdout"),"rU")
        #print "Opened temp file %s for reading."%(os.path.join(tempDir,os.path.basename(procName)+".stdout"))
        stdErrFile = open(os.path.join(tempDir,os.path.basename(procName)+".stderr"),"rU")
        stdOuts[procName] = stdOutFile.read()
        stdErrs[procName] = stdErrFile.read()
        stdOutFile.close()
        stdErrFile.close()
    return stdOuts, stdErrs

def fixFasta(faPath):
    # check a fasta file for requirements, and fix and warn if necessary
    warnedFlag = False
    fa = open(faPath,"r")
    read = fa.read()
    fa.seek(0)
    lines = fa.readlines()
    fa.close()
    if "\r" in read:
        warnedFlag = True
        lines = [line.strip()+"\n" for line in lines]
        log.write("Line endings in fasta file %s converted to unix format\n"%faPath)
    faName = os.path.basename(os.path.splitext(faPath)[0])
    if lines[0].strip()[1:] != faName:
        warnedFlag = True
        lines[0] = ">"+faName+"\n"
        log.write("Replaced title line in fasta file %s with filename\n"%faPath)
    foundUflag = False
    foundCaseFlag = False
    foundUnknownChar = False
    fixedLines = [list(line) for line in lines[1:]]
    for i in xrange(len(fixedLines)):
        for j in xrange(len(fixedLines[i])):
            c = fixedLines[i][j]
            if c=="\n":
                continue
            elif c == "U" or c == "u":
                foundUflag = True
                warnedFlag = True
                fixedLines[i][j] = "T"
            elif c in ["a","t","g","c"]:
                foundCaseFlag = True
                warnedFlag = True
                fixedLines[i][j] = c.upper()
            elif c not in ["A","T","G","C"]:
                foundUnknownChar = True
                warnedFlag = True

    lines = [lines[0]]+["".join(line) for line in fixedLines]
    if foundUflag == True:
        log.write("Replaced U's with T's in fasta file %s\n"%faPath)
    if foundCaseFlag == True:
        log.write("Replaced lowercase sequences with uppercase in fasta file %s\n"%faPath)
    if foundUnknownChar == True:
        log.write("Warning: Found unknown characters (not A,T,G, or C) in fasta file %s.\nThis could cause downstream problems.\n"%faPath)
        
    if warnedFlag == True:
        fa = open(faPath,"w")
        fa.write("".join(lines))
        fa.close()
    
    
#---------------------------------------------------------------------------------
# Parse config files, load config params, create output folders
try:
    log.write("\nStarting config file parsing at %s\n"%(timeStamp()))
    log.flush()

    # parse configuration files                         
    parseConfigFile.parseMain(sys.argv[1])
    # load namespace that now contains configuration parameters
    import conf

    # create output folder if needed
    if not os.path.exists(conf.outputDir):
        os.makedirs(conf.outputDir)

    # check for output sub-folders, create if needed
    folderList = [
          "trimmed_reads",
          "bowtie_index",
          "aligned_reads",
          "counted_mutations",
          "counted_mutations_columns",
          "mutation_strings",
          "mutation_strings_oldstyle",
          "reactivity_profiles",
          "folds"]
    for folder in folderList:
        dirPath = os.path.join(conf.outputDir,folder)
        if not os.path.exists(dirPath):
            os.makedirs(dirPath)

    log.write("Config file successfully parsed at %s\n"%(timeStamp()))
    log.flush()

except Exception as e:
    errorString = "Error:" + str(e) + "\n"
    if "~" not in str(e):  # print a traceback if this is not an error message we created (should really create an Exception subclass)
        errorString = "Error:" + traceback.format_exc() + "\n"
    errorString += "Pipeline initialization failed %s."%timeStamp()
    log.write(errorString+"\n")
    log.flush()
    print errorString
    sys.exit(1)


#-----------------------------------------------------------------------------------------
# Build reference indices as needed
if conf.buildIndex == True:
    try:
        log.write("\nStarting reference index building at %s\n"%(timeStamp()))
        log.flush()

        # link each barcode with a reference index file;
        # avoid duplicating indices
        indexPathDict = {}
        indicesToBuild = {}
        for name in conf.alignments:
            tempPath = "_".join(sorted(conf.alignments[name]["targetNames"]))
            if tempPath not in indexPathDict.values():
                indicesToBuild[tempPath] = conf.alignments[name]["targetNames"]
            indexPathDict[name] = tempPath
        # check fasta sequence titles match filenames, etc. and fix if necessary
        faNames = []
        for subNames in indicesToBuild.values():
            faNames.extend(subNames)
        faNames = set(faNames)
        for faName in faNames:
            faPath = os.path.join(conf.fastaDir,faName+".fa")
            fixFasta(faPath)
            
        # create combination fasta files as needed
        for name in indicesToBuild:
            subNames = indicesToBuild[name]
            if len(subNames) > 1: # single fasta files should already exist, only making combos here
                fastaOut = open(os.path.join(conf.fastaDir,name+".fa"),"w")
                for f in [open(os.path.join(conf.fastaDir,subName+".fa"),"rU") for subName in subNames]:
                    for line in f:
                        fastaOut.write(line)
                    fastaOut.write("\n")
                fastaOut.close()
        
        argListDict = {}
        for filename in indicesToBuild.keys():
            outputPath = os.path.join(conf.outputDir,"bowtie_index/"+filename)
            refPath = os.path.join(conf.fastaDir,filename+".fa")
            argList = ["bowtie2-build", refPath, outputPath]
            argListDict[filename] = argList
        
        maxProc = 4
        stdOuts, stdErrs = spawnProcesses(maxProc, indicesToBuild.keys(), argListDict)

        for filename in stdErrs:
            stdErr = stdErrs[filename]
            stdOut = stdOuts[filename]
            if stdOut != "" and stdOut is not stdOut:
                #print stdOut+"\n"
                pass
            if "Error" in stdErr:
                errorString = "Error: %s bowtie index build failed %s.\n"%(filename,timeStamp())
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("%s bowtie index built successfully %s.\n"%(filename,timeStamp()))
                log.flush()


    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Bowtie index build failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
    


#-------------------------------------------------
# Get list of .fastq sequence files
fileList = os.listdir(conf.inputDir)
fastqFileList = []
for sampleName in conf.alignments.keys():
    fullFileNames = conf.alignments[sampleName]["fileNames"]
    if len(fullFileNames) > 0:
        fastqFileList.extend(fullFileNames)
    else:
        fastqFileList.extend([filename for filename in fileList if os.path.splitext(filename)[1]==".fastq" and
                              "_".join(filename.split("_")[:-4]) == sampleName])
fastqFilePathList = [os.path.join(conf.inputDir,filename) for filename in fastqFileList]


# -------------------------------------------------------------------------------------------
# Run stage 1 (simple read quality filtering)
if conf.trimReads == True:
    try:
        log.write("\nStarting read quality filtering at %s\n"%(timeStamp()))
        log.flush()

        argListDict = {}
        #print "fastqFilePathList: "+str(fastqFilePathList)
        for filePath in fastqFilePathList:
            fileName = os.path.split(filePath)[1]
            argList = [conf.cQualFilterPath, 
                       "-minphred", str(conf.minPhred),
                       "-windowsize", str(conf.windowSize),
                       "-minlength",str(conf.minLength), 
                       "-maxlinelength", "5000",
                       "-filein", filePath,
                       "-fileout", os.path.join(conf.outputDir,"trimmed_reads",fileName)]
            #print " ".join(argList)
            argListDict[filePath] = argList
        maxProc = 4
        stdOuts, stdErrs = spawnProcesses(maxProc, fastqFilePathList, argListDict)

        for filePath in stdErrs:
            fileName = os.path.split(filePath)[1]
            stdErr = stdErrs[filePath]
            stdOut = stdOuts[filePath]
            #if stdOut != "" and stdOut is not None:
            #    print stdOut+"\n"
            #if "Error" in stdErr:
            if len(stdErr) > 0: # TODO: capture exit codes
                errorString = "Error: %s quality filtering failed %s.\n"%(fileName,timeStamp())
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("%s quality filtered successfully.\n"%(fileName))
                log.flush()

    
    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Quality filtering failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
                         
    
#------------------------------------------------------------------------------------------------
# Run stage 2 (sequence alignment using bowtie2)

# By Illumina filenaming convention, fields are underscore delimited.
if conf.alignReads == True:
    samFilePathList = [] # this list will hold SAM file paths so we can sort them in the step following alignment
    try:
        log.write("\nStarting sequence alignment at %s\n"%(timeStamp()))
        log.flush()

        sampleDict = {}
        for sampleName in conf.alignments:
            sampleDict[sampleName] = {"R1":[],"R2":[]}
            fileNames = conf.alignments[sampleName]["fileNames"]
            if len(fileNames) > 0:
            # one or two full filenames were specified in the config file
                sampleDict[sampleName]["R1"] = [os.path.join(conf.outputDir,"trimmed_reads/",fileNames[0])]
                if len(fileNames) > 1:
                    sampleDict[sampleName]["R2"] = [os.path.join(conf.outputDir,"trimmed_reads/",fileNames[1])]
            else:
            # explicit filenames were not specified - find files matching sample name
                for filename in fastqFileList:
                    if filename.split("_")[0] == sampleName: # this effectively disallows underscores in sample names
                    # found a file with a name matching the given sample name
                        readNumberField = filename.split("_")[-2]
                        if readNumberField not in ["R1","R2"]:
                            raise Exception("Error:~ \"%s\" does not match expected filename format.\n"%filename)
                        sampleDict[sampleName][readNumberField] = [os.path.join(conf.outputDir,"trimmed_reads/",filename)]         
        for sampleName in sampleDict.keys():
            sampleDict[sampleName]["R1"].sort()
            sampleDict[sampleName]["R2"].sort()


        indexDir = os.path.join(conf.outputDir,"bowtie_index/")
        argListDict = {}
        for sampleName in sampleDict.keys():
            # exit with warning if space(s) in file path - bowtie2 wrapper can't handle them, at least unescaped
            if " " in sampleDict[sampleName]["R1"]:
                raise Exception("~Error: At least one space found in directory path. Bowtie2/this pipeline can't handle these. Try replacing spaces with underscores in folder names.")
            refPath = os.path.join(indexDir,
                                   "_".join(sorted(conf.alignments[sampleName]["targetNames"])))
            argList = ["bowtie2"]
            
            argList.extend("--local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1".split(" "))

            # param to allow for deletions up to ~200 nucs            
            argList.extend(["--dpad","100"]) 
            # param to allow for deletions up to ~2000 nucs
            #argList.extend(["--dpad","1000"])
            
            argList.extend(["--maxins",str(conf.maxInsertSize)])
            
            argList.extend(["-p","4"]) # use 4 processes for each bowtie instance
            argList.extend(["-x",refPath])
            if conf.alignPaired == True:
                argList.extend(["-1",",".join(sampleDict[sampleName]["R1"])])
                argList.extend(["-2",",".join(sampleDict[sampleName]["R2"])])
            else:
                argList.extend(["-U",",".join(sampleDict[sampleName]["R1"]+sampleDict[sampleName]["R2"])])

            #argList.append("--reorder") # this option seems to cause bowtie to only output R1 reads when also using "-U"
            # I sort the files afterward anyway, so this shouldn't matter
            samFileName = os.path.join(conf.outputDir,"aligned_reads/",sampleName+".sam")
            argList.extend(["-S",samFileName])
            samFilePathList.append(samFileName)
            
            argListDict[sampleName] = argList
            
            #print " ".join(argListDict[sampleName])+"\n"

        maxProc = 1 # align files serially - This does not take full advantage of the parallel capabilities of a computing
                    # cluster. I decided to keep the pipeline simple to run on a desktop computer for the time being.
        # bowtie runs from a wrapper script (so an additional process will be created)
        stdOuts, stdErrs = spawnProcesses(maxProc, sampleDict.keys(), argListDict)
            
        # note - bowtie2 errors output to stderr; alignment stats also output to stderr, so handle this carefully. This may break with future versions of bowtie2
        for sampleName in stdErrs:
            stdErr = stdErrs[sampleName]
            if "Error" in stdErr:
                errorString = "Error: sample %s read alignment failed.\n"%(sampleName)
                errorString += "\n".join(stdErr.splitlines()) # used to only output the last line [-1], but this doesn't catch all errors
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("Sample %s aligned successfully.\n"%(sampleName))
                log.write("Sample %s alignment stats:\n"%sampleName)
                for line in stdErr.splitlines():
                    if "Warning" not in line:
                        log.write("    "+line.strip()+"\n")
                log.flush()

    except Exception as e:
        errorString = str(e)
        # print traceback if this is not a custom exception (should really define an Exception subclass)
        if "~" not in str(e):
            errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Read alignment failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
    
    # sort .sam alignment files so mate-pairs are adjacent
    #try:
    #    log.write("\nStarting alignment file sorting at %s\n"%(timeStamp()))
    #    log.flush()
    
    #    if conf.alignPaired == True:
    #        argListDict = {}
    #        for samFilePath in samFilePathList:
    #            tempFilePath = samFilePath+".sort"
    #            argList = ["sort","-s","-k1,1"] # stable sort on the first key (cluster name)
    #            argList.extend([samFilePath,"-o",samFilePath])
    #            argListDict[samFilePath] = argList

    #        maxProc = 1 # run serially to avoid maxing out node memory in case files are large
    #        stdOuts, stdErrs = spawnProcesses(maxProc, samFilePathList, argListDict)
            
    #        for samFilePath in stdErrs:
    #            stdErr = stdErrs[samFilePath]
    #            if "Error" in stdErr:
    #                errorString = "Error: %s alignment file sorting failed.\n"%(samFilePath)
    #                errorString += "\n".join(stdErr.splitlines()[:1])
    #                log.write(errorString+"\n")
    #                log.flush()
    #                print errorString
    #                sys.exit(1)
    #            else:
    #                log.write("%s alignment file sorted successfully.\n"%(samFilePath))
    #                log.flush()
    #except Exception as e:
    #    errorString = str(e)
    #    # print traceback if this is not a custom exception (should really define an Exception subclass)
    #    if "~" not in str(e):
    #        errorString = "Error:" + traceback.format_exc() + "\n"
    #    errorString += "Alignment file sorting failed %s."%timeStamp()
    #    log.write(errorString+"\n")
    #    log.flush()
    #    print errorString
    #    sys.exit(1)

#---------------------------------------------------------------------------------------------
# Run stage 3 (parse .sam files)

# TODO: check for Seg Faults in subprocesses (they are silently ignored now)
if conf.parseAlignments == True:
#if False:
    try:
        log.write("\nStarting alignment parsing at %s\n"%(timeStamp()))
        log.flush()

        # get list of .sam alignment files
        alignmentDir = os.path.join(conf.outputDir,"aligned_reads/")
        alignedFileList = os.listdir(alignmentDir)
        sampleNameList = [os.path.splitext(filename)[0] \
                          for filename in alignedFileList \
                          if os.path.splitext(filename)[1]==".sam" \
                          and os.path.splitext(filename)[0] in conf.alignments.keys()]
        samFilePathDict = {} # sample name:aligned file path
        indexPathDict = {} # sample name:reference fasta file path
        for sampleName in sampleNameList:
            samFilePathDict[sampleName] = os.path.join(alignmentDir,sampleName+".sam")
            indexPathDict[sampleName] = os.path.join(conf.fastaDir,
                                                     "_".join(sorted(conf.alignments[sampleName]["targetNames"]))+".fa")

        argListDict = {}
        for sampleName in sampleNameList:
            filePath = samFilePathDict[sampleName]
            options = ["-combine_strands","-deletion_masking","-randomly_primed","-trim_both_ends"]
            selected = [conf.alignPaired, conf.longDeletionMasking, conf.randomlyPrimed, conf.trimBothEnds]
            argList = [conf.cParseAlignmentPath]
            for i in range(len(options)):
                if selected[i] == True:
                    argList.append(options[i])
            argList += ["-primer_length", str(conf.primerLength),
                        "-min_map_qual", str(conf.minMapQual),
                        "-ref_seqs", indexPathDict[sampleName],
                        "-file_in", filePath,
                        "-out_folder", os.path.join(conf.outputDir,"mutation_strings/")]
            if conf.removeAmbigDel==True:
                argList += ["-remove_ambig_del"]
            argListDict[sampleName] = argList

        #for sampleName in argListDict.keys():
        #    print "argListDict[sampleName] : "+" ".join(argListDict[sampleName])
            
        maxProc = 6
        stdOuts, stdErrs = spawnProcesses(maxProc, sampleNameList, argListDict)

        for sampleName in stdErrs:
            stdErr = stdErrs[sampleName]
            stdOut = stdOuts[sampleName]
            if stdOut != "" and stdOut is not None:
                print str(stdOut)+"\n"
            if len(stdErr) > 0:
                errorString = "Error: sample %s alignment parsing failed.\n"%(sampleName)
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("Sample %s alignment parsed successfully.\n"%(sampleName))
                log.flush()

    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Alignment parsing failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
        
#-------------------------------------------------------------------------------------------
# Run stage 4 (count mutations)
if conf.countMutations == True:
    try:

        log.write("\nStarting mutation counting at %s\n"%(timeStamp()))
        log.flush()

        # get list of parsed mutation string files
        mutationStringDir = os.path.join(conf.outputDir,"mutation_strings/")
        folderFileList = os.listdir(mutationStringDir)

 
        argListDict = {}
        txtFileList = []
        for sampleName in conf.alignments.keys():
            for target in conf.alignments[sampleName]["targetNames"]:
                fileName = sampleName+"_"+target+".txt"
                if fileName in folderFileList:
                    filePath = os.path.join(mutationStringDir,fileName)
                    txtFileList.append(fileName)
                    refPath =  os.path.join(conf.fastaDir,
                                "_".join(sorted(conf.alignments[sampleName]["targetNames"]))+".fa")
                    outPath = os.path.join(conf.outputDir,"counted_mutations",sampleName+"_"+target+".csv")
        
                    argList = [conf.cCountMutationsPath,
                               "-file_in", filePath,
                               "-ref_seqs", refPath,
                               "-sample_name", sampleName,
                               "-target_name", target,
                               "-file_out", outPath,
                               "-min_phred", str(conf.minPhredToCount)]
                    argListDict[fileName] = argList

        #for fileName in argListDict.keys():
        #    print "argListDict[sampleName] : "+" ".join(argListDict[fileName])
            
        maxProc = 6
        stdOuts, stdErrs = spawnProcesses(maxProc, txtFileList, argListDict)

        for fileName in stdErrs:
            stdErr = stdErrs[fileName]
            stdOut = stdOuts[fileName]
            if stdOut != "" and stdOut is not None:
                print str(stdOut)+"\n"
            if len(stdErr) > 0:
                errorString = "Error: file %s mutation counting failed.\n"%(fileName)
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("File %s mutations counted successfully.\n"%(fileName))
                log.flush()

    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Mutation counting failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)

#---------------------------------------------------------------------------------
# Make mutation strings compatible with older downstream scripts (RING-MaP, not
# included in this distribution)
# Also filter out mutations not meeting phred score cutoff
if conf.countMutations == True and conf.makeOldMutationStrings == True:
    try:
        log.write("\nStarting making old-style mutation strings at %s\n"%(timeStamp()))
        log.flush()

        # get list of parsed mutation string files
        mutationStringDir = os.path.join(conf.outputDir,"mutation_strings/")
        folderFileList = os.listdir(mutationStringDir)

        argListDict = {}
        txtFileList = []
        for sampleName in conf.alignments.keys():
            for target in conf.alignments[sampleName]["targetNames"]:
                fileName = sampleName+"_"+target+".txt"
                if fileName in folderFileList:
                    filePath = os.path.join(mutationStringDir,fileName)
                    txtFileList.append(fileName)
                    refPath =  os.path.join(conf.fastaDir,
                                "_".join(sorted(conf.alignments[sampleName]["targetNames"]))+".fa")
                    outPath = os.path.join(conf.outputDir,"mutation_strings_oldstyle",sampleName+"_"+target+".txt")
        
                    argList = ["python",
                               conf.makeOldMutationStringsPath,
                               filePath,
                               refPath,
                               outPath,
                               str(conf.minPhredToCount)]
                    argListDict[fileName] = argList

        #for fileName in argListDict.keys():
        #    print "argListDict[sampleName] : "+" ".join(argListDict[fileName])
            
        maxProc = 4
        stdOuts, stdErrs = spawnProcesses(maxProc, txtFileList, argListDict)

        for fileName in stdErrs:
            stdErr = stdErrs[fileName]
            stdOut = stdOuts[fileName]
            if stdOut != "" and stdOut is not None:
                print str(stdOut)+"\n"
            if len(stdErr) > 0:
                errorString = "Error: file %s making old mutation strings failed.\n"%(fileName)
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("File %s made old mutation strings successfully.\n"%(fileName))
                log.flush()

    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Making old mutation string failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)

#----------------------------------------------------------------------------------------
# Swap rows to columns for mutation count .csv files (kludgy but useful)
if conf.pivotCSVs == True:
    log.write("\nCopying counted mutation csv's rows to columns at %s\n"%(timeStamp()))
    log.flush()

    # TODO: error handling for this module
    countedDir = os.path.join(conf.outputDir,"counted_mutations")
    colDir = os.path.join(conf.outputDir,"counted_mutations_columns")
    fileList = os.listdir(countedDir)
    countedFileList = [filename for filename in fileList if os.path.splitext(filename)[1]==".csv"]
    countedFilePathList = [os.path.join(countedDir,filename) for filename in countedFileList]
    for csvFilePath in countedFilePathList:
        pivotCSV.pivot(csvFilePath,colDir)

#---------------------------------------------------------------------------------------
# Generate final reactivity profiles

if conf.makeProfiles == True:
    log.write("\nStarting reactivity profile creation at %s\n"%(timeStamp()))
    log.flush()

    # get list of .csv mutation count files
    csvDir = os.path.join(conf.outputDir,"counted_mutations/")

    try:
        argListDict = {}
        for profileName in conf.profiles.keys():
            argList = ["python",
                       conf.genProfilesPath,
                       csvDir,
                       profileName,
                       conf.profiles[profileName]["target"]]
            argList += ["--plus"]
            argList += conf.profiles[profileName]["plus_reagent"]
            argList += ["--minus"]
            argList += conf.profiles[profileName]["minus_reagent"]
            argList += ["--denat"]
            argList += conf.profiles[profileName]["denat_control"]
            #print "profile %s args: %s"%(profileName,str(argList))
            argListDict[profileName] = argList

        maxProc = 4
        stdOuts, stdErrs = spawnProcesses( maxProc, conf.profiles.keys(), argListDict)

        for profileName in stdErrs:
            stdErr = stdErrs[profileName]
            stdOut = stdOuts[profileName]
            if stdOut != "" and stdOut is not None:
                print stdOut+"\n"
            if "Error" in stdErr:
                errorString = "Error: reactivity profile %s generation failed.\n"%(profileName)
                errorString += "\n".join(stdErr.splitlines())
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("Reactivity profile %s generated successfully.\n"%(profileName))
                log.flush()

    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Reactivity profile creation failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
        



#------------------------------------------------------------------------------------------
# Fold RNAs using RNAstructure
if conf.foldSeqs == True:

    try:
        log.write("\nStarting folding at %s\n"%(timeStamp()))
        log.flush()

        seqDir = os.path.join(conf.outputDir,"folds")
        shapeDir = os.path.join(conf.outputDir,"reactivity_profiles")

        # get list of .shape reactivity files
        fileList = os.listdir(shapeDir)
        shapeFileList = [filename for filename in fileList if os.path.splitext(filename)[1]==".shape"]
        
        # create seq files for RNASructure from fasta files
        seqList = []
        for profileName in conf.folds:
            seqName = conf.folds[profileName]
            if seqName not in seqList:
                seqList.append(seqName)
                seqOut = open(os.path.join(seqDir,seqName+".seq"),"w")
                f = open(os.path.join(conf.fastaDir,seqName+".fa"),"rU")
                f.readline()
                seqOut.write(";required comment line\n")
                seqOut.write(seqName.replace(";","")+"\n")
                for line in f:
                    seqOut.write(line.strip().replace("T","U"))
                seqOut.write("1\n")
                seqOut.close()    


        argListDict = {}
        for foldName in conf.folds.keys():        
            seqName = conf.folds[foldName]
            seqPath = os.path.join(seqDir,seqName+".seq")
            ctPath = os.path.join(seqDir,foldName+".ct")
            shapePath =os.path.join(shapeDir,foldName+".shape")
            argList = ["Fold"]
            argList.append(seqPath)
            argList.append(ctPath)
            argList+=["-sh",shapePath]
            argList+=["-sm",str(conf.slope)]
            argList+=["-si",str(conf.intercept)]
            #print str(argList)
            argListDict[foldName] = argList

        maxProc = 4
        stdOuts, stdErrs = spawnProcesses(maxProc, conf.folds.keys(), argListDict)


        for foldName in stdErrs:
            stdErr = stdErrs[foldName]
            stdOut = stdOuts[foldName]
            #if stdOut != "":
            #    print stdOut
            if stdErr != "" and stdErr is not None:
                errorString = "Error: %s folding failed.\n"%(foldName)
                errorString += stdErr
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("%s folded successfully.\n"%(foldName))
                log.flush()
            
    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Folding failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
        
#----------------------------------------------------------------------------------
# Render folded structures using pvclient and PseudoViewer
if conf.renderStructures == True:

    try:
        log.write("\nStarting structure rendering at %s\n"%(timeStamp()))
        log.flush()

        seqDir = os.path.join(conf.outputDir,"folds")
        shapeDir = os.path.join(conf.outputDir,"reactivity_profiles")

        argListDict = {}
        for foldName in conf.folds.keys():       
            ctPath = os.path.join(seqDir,foldName+".ct")
            shapePath = os.path.join(shapeDir,foldName+".shape")
            outPath = os.path.join(seqDir,foldName)
            argList = [conf.renderPath] # pvclient.py in pipeline dir
            argList += ["--ct",ctPath]
            argList += ["--shape",shapePath]
            argList += ["--out",outPath]
            argList += ["--title",foldName]
            #print "argList: "+str(argList)
            argListDict[foldName] = argList

        maxProc = 1 # run serially so we don't overload the Pseudoviewer server
        stdOuts, stdErrs = spawnProcesses(maxProc, conf.folds.keys(), argListDict)

        for foldName in stdErrs:
            stdErr = stdErrs[foldName]
            stdOut = stdOuts[foldName]
            #if stdOut != "":
            #    print stdOut
            if "Error" in stdErr:
                errorString = "Error: %s structure rendering failed.\n"%(foldName)
                errorString += stdErr
                log.write(errorString+"\n")
                log.flush()
                print errorString
                sys.exit(1)
            else:
                log.write("%s structures rendered successfully.\n"%(foldName))
                log.flush()

    except:
        errorString = "Error:" + traceback.format_exc() + "\n"
        errorString += "Structure rendering failed %s."%timeStamp()
        log.write(errorString+"\n")
        log.flush()
        print errorString
        sys.exit(1)
        
doneMessage = "Pipeline completed %s."%timeStamp()
log.write(doneMessage+"\n")
log.flush()
print doneMessage
