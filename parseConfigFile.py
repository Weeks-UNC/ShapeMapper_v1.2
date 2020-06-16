# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Module for parsing configuration file options and storing to easily-loaded file.
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

import sys, os, traceback, pickle

trueVals = ["on","ON","On","True","TRUE","T","true"]
falseVals = ["off","OFF","Off","False","FALSE","F","false"]

def removeComments(line):
# remove comments from line
    poundIndex = line.find("#")
    if poundIndex != -1:
        line = line[:poundIndex]
    return line    

def parseSimpleParams(configFilePath):
# given a path to a config file, parse simple parameters into a dictionary and return it
    params = {}
    configFile = open(configFilePath, "rU")
    
    for line in configFile:
        line = removeComments(line)
        # ignore anything beyond the initial parameter section of the config file
        if line.strip() == "[alignments]":
            break
        splitLine = [c.strip() for c in line.split("=")]
        if len(splitLine) < 2:
            continue
        # if quoted, store string
        if splitLine[1][0] in ["\"","\'"]:
            params[splitLine[0]] = splitLine[1][1:-1]
        # if obvious boolean, store value
        elif splitLine[1] in trueVals:
            params[splitLine[0]] = True
        elif splitLine[1] in falseVals:
            params[splitLine[0]] = False
        # otherwise, convert to number and store
        else:
            if "." in splitLine[1]:
                try:
                    params[splitLine[0]] = float(splitLine[1])
                except ValueError:
                    params[splitLine[0]] = "error"
            else:
                try:
                    params[splitLine[0]] = int(splitLine[1])
                except ValueError:
                    params[splitLine[0]] = "error"
    configFile.close()
    return params

def parseOtherParams(configFilePath):
# given a file handle, parse other options (such as alignment targets)
    alignments = {}
    profiles = {}
    folds = {}
    section = ""
    profileName = ""
    configFile = open(configFilePath, "rU")
    for line in configFile:
        line = removeComments(line)
        stripLine = line.strip()
        if stripLine in ["[alignments]","[profiles]","[folds]"]:
            section = stripLine[1:-1]
            continue
        if section == "alignments":
            splitLine = [s.strip() for s in stripLine.split("=")]
            if len(splitLine) > 1:
                leftSplit = [s.strip() for s in splitLine[0].split(":")]
                targetNames = [s.strip() for s in splitLine[1].split(",")]
                sampleName = leftSplit[0]
                fileNames = []
                if len(leftSplit) > 1:
                    fileNames = [s.strip() for s in leftSplit[1].split(",")]
                alignments[sampleName] = {"fileNames":fileNames,"targetNames":targetNames}
        elif section == "profiles":
            splitLine = [s.strip() for s in stripLine.split("=")]
            if len(splitLine) > 1:
                varName = splitLine[0]
                value = splitLine[1]
                if varName == "name":
                    profileName = value
                    profiles[profileName] = {}
                elif varName == "target":
                    profiles[profileName]["target"] = value
                elif varName in ["plus_reagent","minus_reagent","denat_control"]:
                    valueList = [value.strip() for value in value.split(",")]
                    profiles[profileName][varName] = valueList
        elif section == "folds":
            foldName = line.strip()
            if foldName != "":
                folds[foldName] = ""
    #print "alignments:"+str(alignments)
    #print "profiles:"+str(profiles)
    #print "folds:"+str(folds)
    for foldName in folds.keys():
        if foldName in profiles.keys():
            folds[foldName] = profiles[foldName]["target"]
        else:
            folds[foldName] = "no matching profile found"
    return alignments, profiles, folds



def checkSimpleParams(defaultParams, userParams):
# compare user-supplied parameters to defaults
# generate error messages if there are apparent typos
    for key in userParams.keys():
        if key not in defaultParams.keys():
            raise Exception("~ Error in %s: Unrecognized parameter \"%s\""%(sys.argv[1],key))
        userType = type(userParams[key])
        defaultType = type(defaultParams[key])
        if userParams[key] == "error":
            raise Exception("~ Error in %s: Parameter \"%s\" incorrectly formatted."%(sys.argv[1],key))
        elif userType != defaultType:
            raise Exception("~ Error in %s: Parameter \"%s\" accepts %s, but got %s."%(sys.argv[1],key,defaultType,userType))

def matchFileNames(name, fileList, requireFullMatch=True):
    foundName = False
    for fileName in fileList:
        if requireFullMatch == True:
            if fileName == name:
                foundName = True
        else:
            if fileName[:len(name)] == name:
                foundName = True
    if foundName == False:
        raise Exception("~ Error in %s: no fastq found matching \"%s\"."%(sys.argv[1], name))
    
def checkOtherParams(alignments, profiles, folds, paths):
    # check if alignment sample names match existing files
    fileList = os.listdir(paths["inputDir"])
    fastqFileList = [filename for filename in fileList if os.path.splitext(filename)[1]==".fastq"]
    for name in alignments.keys():
        f = alignments[name]["fileNames"]
        if len(f) > 0:
            for fullName in f:
                matchFileNames(fullName,fastqFileList)
        else:
            matchFileNames(name,fastqFileList,requireFullMatch=False)
    # check that alignment targets match existing fasta files
    fastaFileList = [filename for filename in fileList if os.path.splitext(filename)[1]==".fa"]
    for name in alignments.keys():
        for target in alignments[name]["targetNames"]:
            foundFasta = False
            for fileName in fastaFileList:
                if os.path.splitext(fileName)[0] == target:
                    foundFasta = True
            if foundFasta == False:
                raise Exception("~ Error in %s: no fasta (.fa) found matching \"%s\"."%(sys.argv[1], target))
    # check if profile dictionary matches basic expectations
    # target name matches alignments
    # sample names are present in alignments
    for name in profiles.keys():
        profile = profiles[name]
        target = profile["target"]
        conditions = [profile["plus_reagent"],profile["minus_reagent"],profile["denat_control"]] 
        for condition in conditions:
            for sampleName in condition:
                foundSampleName = False
                foundTarget = False
                if sampleName in alignments.keys():
                    foundSampleName = True
                    if target in alignments[sampleName]["targetNames"]:
                        foundTarget = True
                if foundSampleName == False:
                    raise Exception("~ Error in %s: reactivity profile name \"%s\": sample name \"%s\" does not match any name in alignments section."%(sys.argv[1], name, sampleName)) 
                if foundTarget == False:
                    raise Exception("~ Error in %s: reactivity profile name \"%s\": target name \"%s\" does not match target listed in alignments section."%(sys.argv[1], name, target))

    # check if names in folds match names in profiles
    for name in folds.keys():
        if name not in profiles.keys():
            raise Exception("~ Error in %s: fold name \"%s\" does not match any reactivity profile name."%(sys.argv[1],name))


def parseMain(userConfigPath):
    #-----------------------------
    # locations of "modules" to spawn
    paths = {}
    thisScriptPath = os.path.realpath(__file__)
    thisScriptDir = os.path.split(thisScriptPath)[0]
    paths["cQualFilterPath"] = os.path.join(thisScriptDir,"trimPhred")
    paths["countMutationsPath"] = os.path.join(thisScriptDir,"countMutations.py")
    paths["cParseAlignmentPath"] = os.path.join(thisScriptDir,"parseAlignment")
    paths["cCountMutationsPath"] = os.path.join(thisScriptDir,"countMutations")
    paths["genProfilesPath"] = os.path.join(thisScriptDir,"generateReactivityProfiles.py")
    paths["renderPath"] = os.path.join(thisScriptDir,"pvclient.py")
    paths["makeOldMutationStringsPath"] = os.path.join(thisScriptDir, "makeOldMutationStrings.py")

    #---------------------------------------------------------
    # default directories for data input, output, and reference sequences
    # TODO: expose config file options for using pre-built reference index (useful for aligning to large genomes)
    # input dir is current working directory
    paths["inputDir"] = os.getcwd() 
    # output dir is "output" subfolder within current working directory 
    paths["outputDir"] = os.path.join(paths["inputDir"],"output")
    # location of reference fasta files defaults to the input directory
    paths["fastaDir"] = paths["inputDir"]

    alignments = {}
    profiles = {}
    folds = {}


    # load default parameters from defaults.cfg
    simpleDefParams = parseSimpleParams(os.path.join(thisScriptDir,"defaults.cfg"))
    #print "simpleDefParams:"+str(simpleDefParams)

    # load user parameters
    simpleUserParams = parseSimpleParams(userConfigPath)

    # check user params for errors
    checkSimpleParams( simpleDefParams, simpleUserParams)

    # load other parameters
    alignments, profiles, folds = parseOtherParams(userConfigPath)

    # check other params for errors
    if simpleUserParams.get("devMode",False) == False:
        checkOtherParams(alignments, profiles, folds, paths)

    # put all config options in a single dictionary
    conf = {}
    for key in simpleDefParams.keys():
        conf[key] = simpleDefParams[key]
    for key in simpleUserParams.keys():
        conf[key] = simpleUserParams[key]
    for key in paths.keys():
        conf[key] = paths[key]
    conf["alignments"] = alignments
    conf["profiles"] = profiles
    conf["folds"] = folds

    # write conf to file - pipeline modules will import a helper script called conf to get access to namespace
    fileOut = open(os.path.join(os.getcwd(),"temp_config.pickle"), "w")
    pickle.dump(conf, fileOut)
    fileOut.close()


