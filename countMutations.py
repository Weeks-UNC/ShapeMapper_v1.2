# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Counts up mutations from sequencing data.
# Copyright Steven Busan 2014

#---------------------------------------------------------------------------------------------------
# GPL statement:
#
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
#
#----------------------------------------------------------------------------------------------------

import sys, os, traceback, copy, re, time

try:
    import conf
    import numpy
    import matplotlib.pyplot as plot

    def safeDivide(num,den):
    # convert input numbers to floats and divide
    # If denominator is 0, output -999 to indicate error.
        result = -999.0
        try:
            result = float(num)/float(den)
        except ZeroDivisionError:
            pass
        return result
    
    def parseFasta(fastaFile):
    # given an open fasta file, return a dict of sequence names
    # and associated sequence strings
        seq = {}
        lines = fastaFile.readlines()
        currentSeqName = ""
        for line in lines:
            if line[0] == ">":
                # sequence name found
                currentSeqName = line.strip()[1:]
                seq[currentSeqName] = ""
            else:
                # remove whitespace and append to current sequence
                seq[currentSeqName] += "".join(line.split())
        return seq


    def parseCigarString(cigarString, read, quals, refSeq, startIndex):
        # define regular expression patterns for parsing CIGAR string
        item = re.compile(r"[0-9]+[M|I|D|S]")
        number = re.compile(r"([0-9]+)")
        splitCigarString = item.findall(cigarString.strip())
        cigarList = [number.split(x)[1:] for x in splitCigarString]
        #print "CIGAR: %s\nsplit CIGAR: %s"%(cigarString,str(cigarList))
        readCopy = read        
        alignedQuals = {}
        lastQual= "#"
        refIndex = startIndex # begin at leftmost matching position
        events = {} # keys will be nuc indices with respect to reference sequence
        #inserts = {}
        for i in range(len(cigarList)):
            region = cigarList[i]
            regionLength = int(region[0])
            regionType = region[1]
            # alignment match (could be match or mismatch, use raw read to resolve)
            if regionType == "M":
                matchRegion = read[:regionLength]
                qualsRegion = quals[:regionLength]
                read = read[regionLength:]
                quals = quals[regionLength:]
                for regionIndex in range(len(matchRegion)):
                    nuc = matchRegion[regionIndex]
                    qual = qualsRegion[regionIndex]
                    lastQual = qual
                    alignedQuals[refIndex] = qual
                    if nuc != refSeq[refIndex]:
                        # nuc mismatch found
                        #print "Mismatch at ref index %i"%refIndex
                        events[refIndex] = nuc
                    else:
                        events[refIndex] = "|"
                    refIndex += 1
            # insertion
            elif regionType == "I":
                insertedSeq = read[:regionLength]
                #inserts[refIndex] = insertedSeq
                read = read[regionLength:]
                quals = quals[regionLength:]
                # ignore for now (insertions are usually unhelpful)
            # deletion
            elif regionType == "D":
                #print "Deletion %i nucs long at ref index %i"%(regionLength, refIndex)
                for deletionIndex in range(regionLength):
                    events[refIndex] = "-"
                    alignedQuals[refIndex] = lastQual # don't have a phred score for a deletion (no basecall associated with position), so just use last nearby score
                    refIndex += 1
            # padding
            elif regionType == "S":
                read = read[regionLength:]
                quals = quals[regionLength:] # missing from v12b and all pipeline versions before - results in misaligned phred scores and errors combining read pairs
                if i == len(cigarList)-1: # rightmost end of read
                    for offsetIndex in range(regionLength):
                        padIndex = refIndex+offsetIndex
                        if padIndex >= 0 and padIndex < len(refSeq):
                        #if True:
                            events[padIndex] = "s"
                            alignedQuals[padIndex] = "#"
                elif i == 0: # leftmost end of read 
                    for offsetIndex in range(regionLength+1):
                        padIndex = refIndex-offsetIndex
                        if padIndex > 0 and padIndex < len(refSeq):
                        #if True:
                            events[padIndex] = "s"
                            alignedQuals[padIndex] = "#"
        sortedKeys = sorted(events.keys())
        printEvents = ""
        printQuals = ""
        for i in xrange(min(sortedKeys),max(sortedKeys)+1):
            try:
                printEvents += events[i]
            except KeyError:
                printEvents += " "
            try:
                printQuals += alignedQuals[i]
            except KeyError:
                printQuals += " "
                    
        print printEvents
        print printQuals
        return events, alignedQuals

    def cullEvents(events, refSeq, seqName, longDeletions):
        eventChars = ["s","-","A","T","G","C","N"]
        culledEvents = dict(events)
        totalMismatches = 0
        totalDeletions = 0
        unambiguousDeletions = 0
        # scan left to right
        #sortedKeys = sorted(events.keys(), reverse=True)
        for i in range(min(events.keys()), max(events.keys())+1):
        #for i in sortedKeys
            if i not in events.keys():
                culledEvents[i] = "~"
            elif events[i] in ["A","T","G","C"] and events[i+1] not in eventChars:
                # mismatch with no adjacent downstream mutation
                culledEvents[i] = events[i]
                totalMismatches += 1
            elif events[i] == "-":
                if events[i+1] not in eventChars:
                # deleted nuc with no adjacent downstream mutations
                    totalDeletions += 1
                    # ensure deletion is 'unambiguously' aligned (this isn't guaranteed in highly repetitive regions,
                    # but this should take care of the vast majority of ambiguous deletions)
                    k = 1
                    while events[i-k] == "-":
                        k += 1
                    fivePrimeStart = i-k+1
                    threePrimeEnd = i
                    ambiguous = False
                    # slide deletion upstream and downstream
                    notDelSeq = refSeq[:fivePrimeStart]+refSeq[threePrimeEnd+1:]
                    maxOffset = 1+threePrimeEnd-fivePrimeStart
                    for offset in range(-maxOffset,maxOffset+1):
                        notSubSeq = ""
                        if offset != 0:
                            try:
                                notSubSeq = refSeq[:fivePrimeStart+offset]+refSeq[threePrimeEnd+offset+1:]
                            except:
                                pass
                        if notSubSeq == notDelSeq:
                            ambiguous = True
                    
                    if ambiguous == False:
                        culledEvents[i] = "-"
                        unambiguousDeletions += 1
                        if conf.outputLongDeletions == True:
                            longDeletions[seqName].setdefault(threePrimeEnd+1,{}).setdefault(fivePrimeStart+1,0)
                            longDeletions[seqName][threePrimeEnd+1][fivePrimeStart+1] += 1
                    else:
                        # remove this deletion
                        culledEvents[i] = "|"
                        # remove any deletion-masked nucleotides if this deletion is flagged as ambiguous
                        try:
                            offset = -1
                            while events[i+offset] == "-":
                                culledEvents[i+offset] = "|"
                                offset -= 1
                        except IndexError:
                            pass
                # this must be a deleted nucleotide not at the 3 prime end of a long deletion
                elif conf.longDeletionMasking == True:
                    # if deletion masking is turned on, this nuc will not contribute to coverage
                    culledEvents[i] = "~"
                else:
                    # if deletion masking is turned off, this nuc will contribute to coverage (somewhat incorrectly, but will be easier to parse for single-molecule analyses)
                    culledEvents[i] = "|"
            elif events[i] not in ["~","|","s"]:
                # if this position has been culled, replace with match character
                culledEvents[i] = "|"
        return [culledEvents, totalMismatches, totalDeletions, unambiguousDeletions]
    
    def parseLine(line):    
        splitLine = line.split()
        # extract useful fields
        clusterName = splitLine[0]
        targetName = splitLine[2]
        if targetName not in refSeqs.keys() or clusterName[0] == "@":
            return [None, None, None, None, None, None]
        # index of leftmost ref nuc in this alignment (converted from 1-based to 0-based index)
        startIndex = int(splitLine[3])-1
        cigarString = splitLine[5]
        mappingPos = int(splitLine[3])
        mappingQual = int(splitLine[4])
        rawRead = splitLine[9]
        rawQuals = splitLine[10]
        if cigarString != "*":
            refSeq = refSeqs[targetName]
            print clusterName
            events, qualities = parseCigarString(cigarString, rawRead, rawQuals, refSeq, startIndex)
            #culledEvents = cullEvents(events, refSeq)
            return [clusterName, targetName, events, qualities, mappingQual, rawRead]
        else:
            return [None, None, None, None, None, None]

    def combineEvents(R1, R2, Q1, Q2):
    # given 2 mutation dicts with nuc indices as keys, combine and return a single dict 
    # also return a list of overlapping nucleotides
        combinedEvents = dict(R1)
        overlappingNucs = []
        for i in R2.keys():
            combinedEvents[i] = R2[i]
            if i in R1.keys():
                overlappingNucs.append(i)
                # handle positions which do not agree between mate pairs
                if R1[i] != R2[i]:
                    try:
                        if R1[i] == "s":
                            combinedEvents[i] = R2[i]
                        elif R2[i] == "s":
                            combinedEvents[i] = R1[i]
                        elif ord(Q2[i])-33 > ord(Q1[i])-33: # favor higher phred scoring nucs
                            combinedEvents[i] = R2[i]
                            # TODO: this could break up long deletions - fix this
                        else: # fall back to read one if quality scores are equal but nucs are different
                            combinedEvents[i] = R1[i]
                    except IndexError:
                        print "len(R1) = %i, len(Q1) = %i, len(R2) = %i, len(Q2) = %i"%(len(R1),len(Q1),len(R2),len(Q2))
                        raise
                else:
                    combinedEvents[i] = R2[i]
        return combinedEvents, overlappingNucs


    filePath = sys.argv[1]
    fileGen = None
    if conf.lowMemoryMode == False:
        fileIn = open(filePath, "rU")
        fileGen = iter(fileIn.readlines())
        fileIn.close()
    else:
        fileGen = iter(open(filePath, "rU"))
    fileName = os.path.split(filePath)[1]
    refPath = sys.argv[2]
    refFile = open(refPath, "rU")
    sampleName = os.path.splitext(fileName)[0]
    fileOutPath = os.path.join(conf.outputDir,"counted_mutations/",sampleName+".csv")
    fileOut = open(fileOutPath, "w")

    # read in reference sequences and names
    refSeqs = parseFasta(refFile)
    refFile.close()

    # create data structures to store event counts and read depths
    nucs = ["A","T","G","C"]
    def rmChar(listIn, charIn):
        return [charOut for charOut in listIn if charOut != charIn]
    mismatch = {}
    deletion = {}
    totalMismatches = {}
    totalDeletions = {}
    unambiguousDeletions = {}
    longDeletions = {}
    insertion = {}
    depth = {}
    doubleCountDepth = {}
    mutationCount = {}
    mutationRate = {}
    primerCounts = {}
    culledEventsFiles = {}
    longDeletionsFiles = {}    
    perReadLengthsUnpaired = {}
    perReadLengthsPaired = {}

    sys.stderr.write("sample %s ref seq names:"%sampleName)
    
    for seqName in refSeqs:
        sys.stderr.write(str(seqName)+"\t")
        #sys.stderr.write(str(len(refSeqs[seqName]))+"\n")
        #print seqName+": "+str(len(refSeqs[seqName]))+"\n"
        mismatch[seqName] = {}
        for fromNuc in nucs:
            mismatch[seqName][fromNuc] = {}
            for toNuc in rmChar(nucs,fromNuc):
                mismatch[seqName][fromNuc][toNuc] = [0]*len(refSeqs[seqName])
        deletion[seqName] = {}
        totalMismatches[seqName] = 0
        totalDeletions[seqName] = 0
        unambiguousDeletions[seqName] = 0
        longDeletions[seqName] = {}
        for fromNuc in nucs:
            deletion[seqName][fromNuc] = [0]*len(refSeqs[seqName])
        insertion[seqName] = {}
        for i in range(len(refSeqs[seqName])):
            insertion[seqName][i] = {"":0}
        depth[seqName] = [0]*len(refSeqs[seqName])
        doubleCountDepth[seqName] = [0]*len(refSeqs[seqName])
        mutationCount[seqName] = [0]*len(refSeqs[seqName])
        mutationRate[seqName] = [0.0]*len(refSeqs[seqName])
        primerCounts[seqName] = [0]*len(refSeqs[seqName])
        if conf.outputMutationStrings == True:
            culledEventsFiles[seqName] = open(os.path.join(conf.outputDir,"mutation_strings/",sampleName+"_"+seqName+".txt"), "w")
        if conf.outputLongDeletions == True:
            longDeletionsFiles[seqName] = open(os.path.join(conf.outputDir,"long_deletions/",sampleName+"_"+seqName+".txt"),"w")
        perReadLengthsUnpaired[seqName] = [0]*500
        perReadLengthsPaired[seqName] = [0]*500
    sys.stderr.write("\n")
    
    covered = ["|","-","A","T","G","C"]
    lineCount = 0
    prevClusterName = ""
    prevTargetName = ""
    prevLine = ""
    prevEvents = []
    prevQualities = []
    prevRawRead = ""
    pairFound = False
    updateCounts = False
    breakLoop = False    

    while breakLoop == False:
        try:
            line = fileGen.next()
            lineCount += 1
            returnData = parseLine(line)
            clusterName = returnData[0]
            targetName = returnData[1]
            events = returnData[2]
            qualities = returnData[3]
            mappingQual = returnData[4]
            rawRead = returnData[5]
        except StopIteration:
            # let the loop run one line past the end of the file to allow the last read to be counted
            # ugly, but it works
            clusterName, targetName, events = None, None, None
            qualities = None
            mappingQual = None
            rawRead = None
            breakLoop = True
        
        isCombinedRead = False
        
        if line[0] == "@":
            # ignore headers
            continue

        # exclude reads whose mapping location is not clearly unique (i.e. this read corresponds to multiple reference sequences)
        if mappingQual < conf.minMapQual:
            continue
                                
        #print line.strip()
        
        if clusterName is not None:                 
            # TODO: replace this awkward control flow with a two-pass system: read through once and index (paired/unpaired), then read through again and parse
            if conf.alignPaired == True:
                #print "clusterName: "+clusterName
                #print "prevCluster: "+prevClusterName
                if pairFound == False:
                    #if clusterName==prevClusterName and targetName != prevTargetName:
                    #    print "    Read target does not match that of the mate pair."
                    if clusterName == prevClusterName and targetName == prevTargetName:
                        print "    Found mate pair with matching target."
                        pairFound = True
                        updateCounts = True
                        isCombinedRead = True
                        prevSortedKeys = sorted(prevEvents.keys())
                        sortedKeys = sorted(events.keys())
                        #print "mate 1 min, max: %i, %i"%(min(prevSortedKeys),max(prevSortedKeys))
                        #print "mate 2 min, max: %i, %i"%(min(sortedKeys),max(sortedKeys))
                        #print "mate 1:"+"".join([prevEvents[ind] for ind in prevSortedKeys])
                        #print "mate 2:"+"".join([events[ind] for ind in sortedKeys])
                        eventsToWrite, overlappingNucs = combineEvents(prevEvents, events, prevQualities, qualities)
                        sortedKeys = sorted(eventsToWrite.keys())
                        #print "combined min, max: %i, %i"%(min(sortedKeys),max(sortedKeys))
                        #print "combined events:"+"".join([eventsToWrite[ind] for ind in sortedKeys])
                        
                        #primerCounts[targetName][max(eventsToWrite.keys())] += 1
                    elif prevClusterName != "":  # need this check, otherwise will try to update counts with blank events from non-existent 0th line
                        if conf.ignoreUnpaired == False:
                            print "    Non-matching read. Need to update counts with previous read info."
                            eventsToWrite = dict(prevEvents)
                            overlappingNucs = []
                            pairFound = False
                            updateCounts = True
                        else:
                            #print "    Non-matching read. Ignoring."
                            pass
                else:
                    #print "    Updated counts last read, so wait one more read."
                    pairFound = False
                    updateCounts = False
            elif prevClusterName != "" and conf.ignoreUnpaired == False:
                updateCounts = True
                eventsToWrite = dict(prevEvents)
                overlappingNucs = []
                sortedKeys = sorted(eventsToWrite.keys())
                #print "".join([eventsToWrite[key] for key in sortedKeys])
                #primerCounts[prevTargetName][max(eventsToWrite.keys())] += 1
            
            if updateCounts == True:
                #print "    Removing primer region."
                # TODO: have random primer clipping option for each target RNA, instead of "global" option   
                # blank out rightmost primer region of read (and 1 further nuc) with "~" char
                if conf.randomlyPrimed==True:
                    endIndex = max(eventsToWrite.keys())
                    for i in range(endIndex, endIndex-conf.primerLength-1, -1):
                        eventsToWrite[i] = "~"
                #print "    Updating counts."
                #print "        prevTargetName: %s"%prevTargetName
                #print "        targetName: %s"%targetName
                updateCounts = False
                refSeq = refSeqs[prevTargetName]
                returned = cullEvents(eventsToWrite, refSeq, prevTargetName, longDeletions) # longDeletions is passed by reference, and will be modified in place
                ROI = range(1,len(refSeq)+1)
                culledEventsToWrite = returned[0]
                mismatchCount = returned[1]
                delCount = returned[2]
                unambigDelCount = returned[3]
                        
                totalMismatches[prevTargetName] += mismatchCount
                totalDeletions[prevTargetName] += delCount
                unambiguousDeletions[prevTargetName] += unambigDelCount
                # write mutation strings to file
                if conf.outputMutationStrings == True:
                    sortedKeys = sorted(culledEventsToWrite.keys())
                    #sortedKeys = [index for index in sortedKeys if index>=0 and index<len(refSeq)] 
                    #print "culled min, max: %i, %i"%(min(sortedKeys),max(sortedKeys))
                    #print "culled events:"+"".join([culledEventsToWrite[ind] for ind in sortedKeys])
                        
                    # format:
                    # <start nuc (1-based)> TAB <end nuc> TAB <events>
                    culledEventsFiles[prevTargetName].write(str(sortedKeys[0]+1)+"\t"+str(sortedKeys[-1]+1)+"\t"+"".join([culledEventsToWrite[key] for key in sortedKeys])+"\n")
                
                #print "".join([culledEventsToWrite[key] for key in sortedKeys])+"\n"
                for i in culledEventsToWrite.keys():
                    c = culledEventsToWrite[i]
                    if c in covered:
                        #try:
                        depth[prevTargetName][i] += 1
                        doubleCountDepth[prevTargetName][i] += 1
                        if i in overlappingNucs:
                            doubleCountDepth[prevTargetName][i] += 1
                        #except IndexError:
                            #print "i %i out of range within depth[%s]"%(i,prevTargetName)
                            #sys.stdout.flush()
                            #time.sleep(20)
                    if c in nucs or c=="-":
                        mutationCount[prevTargetName][i] += 1
                    if c in nucs:
                        #try:
                        mismatch[prevTargetName][refSeq[i]][c][i] += 1
                        #except KeyError:
                        #    print "mismatch key error"
                        #    print "refSeq: %s"%refSeq
                        #    time.sleep(10)
                    elif c == "-":
                        deletion[prevTargetName][refSeq[i]][i] += 1
                    #elif c == "<":
                        #primerCounts[targetName][i] += 1
            prevClusterName = clusterName
            prevTargetName = targetName
            prevLine = line
            prevEvents = events
            prevRawRead = rawRead
            prevQualities = qualities

    # close culled events files
    if conf.outputMutationStrings == True:
        for seqName in refSeqs:
            culledEventsFiles[seqName].close()
 
    # calculate per-nuc mutation frequencies (deletions+mismatches)/depth
    for seqName in refSeqs:
        for i in range(len(depth[seqName])):
            mutationRate[seqName][i] = safeDivide(mutationCount[seqName][i],depth[seqName][i])
        # handle excluded nucs
        for i in range(len(mutationRate[seqName])):
            if depth[seqName][i] == 0:
                mutationRate[seqName][i] = -999.0
    
    
    for seqName in refSeqs:

        # write mutation and coverage counts
        fileOut.write("%s,%s,sequence,%s\n"%(sampleName,seqName,",".join([c for c in refSeqs[seqName]])))
        for toNuc in deletion[seqName]:
            fileOut.write("%s,%s,%s del,%s\n"%(sampleName,seqName,toNuc,",".join([str(n) for n in deletion[seqName][toNuc]])))
        for fromNuc in mismatch[seqName]:
            for toNuc in mismatch[seqName][fromNuc]:
                fileOut.write("%s,%s,%s->%s,%s\n"%(sampleName,seqName,fromNuc,toNuc,",".join([str(n) for n in mismatch[seqName][fromNuc][toNuc]])))
        fileOut.write("%s,%s,mutation rate,%s\n"%(sampleName,seqName,",".join([str(n) for n in mutationRate[seqName]])))
        fileOut.write("%s,%s,depth,%s\n"%(sampleName,seqName,",".join([str(int(n)) for n in depth[seqName]])))
        fileOut.write("%s,%s,depth double-counted paired,%s\n"%(sampleName,seqName,",".join([str(int(n)) for n in doubleCountDepth[seqName]])))
        fileOut.write("%s,%s,total mismatches,%i\n"%(sampleName,seqName,totalMismatches[seqName]))
        fileOut.write("%s,%s,total deletions,%i\n"%(sampleName,seqName,totalDeletions[seqName]))
        fileOut.write("%s,%s,unambiguous deletions,%i\n"%(sampleName,seqName,unambiguousDeletions[seqName]))
        fileOut.write("\n")
    
except Exception as e:
    sys.stderr.write(traceback.format_exc())
    exit(1)

