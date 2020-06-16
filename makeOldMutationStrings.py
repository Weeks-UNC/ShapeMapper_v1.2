# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Convert parsed alignment files (mutation strings) into
# a simpler format for reverse compatibility with some old
# scripts (e.g. RING-MaP, not included in this distribution)
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

import sys, conf

fileIn = open(sys.argv[1], "rU")
fileOut = open(sys.argv[3], "w")

refSeq = "".join([x.strip() for x in open(sys.argv[2], "rU").readlines()[1:]])
seqLen = len(refSeq)

for line in fileIn:
    splitLine = line.strip().split()
    start = int(splitLine[1])
    end = int(splitLine[2])
    read = list(splitLine[3])
    quals = list(splitLine[5])
    if start < 1:
        read = read[-start+1:]
        quals = quals[-start+1:]
        start = 1
    if end > seqLen:
        read = read[:seqLen-end]
        quals = quals[:seqLen-end]
        end = seqLen
    for i in xrange(len(quals)-1):
        mutChars = ['A','T','G','C','-']
        if read[i] in mutChars and read[i+1] not in mutChars:
            # rightmost nuc in a mutation
            foundBadNuc = False
            offset = 0
            while (i-offset>=0) and (read[i-offset] in mutChars):
                if (ord(quals[i-offset])-33) < conf.minPhredToCount:
                    foundBadNuc = True
                offset += 1
            if foundBadNuc == True:
                offset = 0
                while (i-offset>=0) and (read[i-offset] in mutChars):
                    read[i-offset] = '|' #'~'
                    offset += 1
        if read[i] == '~':
            read[i] = '|'
        if quals[i] == '!':
            read[i] = '~'

    fileOut.write("%i\t%i\t%s\n"%(start,end,"".join(read)))
