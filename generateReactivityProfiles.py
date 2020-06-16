# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Make reactivity profiles from mutation count data.
# Copyright Steven Busan and Greggory Rice 2014

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

import sys, os, traceback, copy, re, math

import numpy

rxColor = "red"
bgColor = "blue"
dcColor = "green"

def findArg(search, argList):
    # return first index of arg in arglist matching search string
    for i in range(len(argList)):
        if search == argList[i]:
            return i

def combineData(csvFileList, csvFilePathList, sampleNames):
    # make combined reagent mutation frequencies, coverage depths
    depths = {}
    counts = {}
    rates = {}
    seq = {}
    for n in xrange(len(csvFileList)):
        csvFileName = csvFileList[n]
        csvFilePath = csvFilePathList[n]
        for sampleName in sampleNames:
            if csvFileName.startswith(sampleName)==False:
                continue
            fileIn = open(csvFilePath, "rU")
            print "loading data from file %s"%csvFilePath
            sys.stdout.flush()
            for line in fileIn:
                splitLine = line.strip().split(",")
                if len(splitLine) > 5:
                    lineSampleName = splitLine[0]
                    # Check that the sample name within the file is correct.
                    # Using the filename alone could result in double-counting
                    # in some rare situations if sample names overlap.
                    if lineSampleName not in sampleNames:
                        continue
                    lineTargetName = splitLine[1]
                    lineType = splitLine[2]
                    if lineTargetName == targetName:
                        dataList = splitLine[3:]
                        if "sequence" in lineType:
                            i = 0
                            for item in dataList:
                                if item!='':
                                    seq[i] = item.strip()
                                    i += 1
                                else:
                                    break
                        elif "->" in lineType or (conf.ignoreDeletions==False and "del" in lineType and "deletion" not in lineType):
                            i = 0
                            for item in dataList:
                                if item!='':
                                    counts[i] = counts.get(i,0) + int(item)
                                    i += 1
                                else:
                                    break
                        elif "depth" in lineType and "double" not in lineType:
                            i = 0
                            for item in dataList:
                                if item!='':
                                    depths[i] = depths.get(i,0) + int(item)
                                    i += 1
                                else:
                                    break
            fileIn.close()
            del fileIn
    
    for i in counts.keys():
        try:
            rates[i] = counts[i]/float(depths[i])
        except ZeroDivisionError:
            rates[i] = -999.0

    return rates.values(), depths.values(), seq.values()

#Following 3 functions modified from Gregg Rice's boxplot normalization script
#
# 1.find the scaling factor by ranking  all of the shape values
# 2. take 1.5*abs(Q1-Q3) as a cutoff
# 3. remove either the top 10% of the RNA or the positions above this cutoff, whichever is smaller
# 4. Average the next 10% from the original length of the RNA --> this is the scaling factor
def calcQuartile(x,q,qtype=7):
    #source: http://adorio-research.org/wordpress/?p=125
    # x = array, q = quartile (in % as a decimal)
    y=x
    n = len(y)
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
      (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
      (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
      (0,   0, 0, 1), # California linear interpolation, R type 4
      (0.5, 0, 0, 1), # hydrologists method, R type 5
      (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
      (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
      (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
      (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
     ]
    a, b, c, d = abcd[qtype-1]
    g, j = math.modf( a + (n+b) * q -1)
    if j<0:
        return x[0]
    elif j>=n:
        return x[n-1]
    j = int(math.floor(j))
    if g == 0:
        return x[j]
    else:
        return y[j] + (y[j+1]- y[j])* (c + d * g)

def findBoxplotFactor(array):
    x,o,a = [],[],0
    #x = [n for n in array if n > -500]
    # Following line is behavior that normalization and structure modeling were optimized with,
    # although this behavior is probably not ideal
    x = [n if n>-500 else 0 for n in array]
    if len(x)/10 < 1:
        normFactor = 1.0
    else:
        x.sort()
        tenPct = len(x)/10
        fivePct = len(x)/20
        #calculate the interquartile range *1.5
        qLimit = 1.5*abs(calcQuartile(x,0.25)-calcQuartile(x,0.75))
        tenLimit = x[len(x)-1 - tenPct]
        fiveLimit = x[len(x)-1 - fivePct]
        #choose the cutoff that eliminates the fewest points
        limit = max(qLimit,tenLimit)
        if len(x)<100:
            limit = max(qLimit,fiveLimit)
        #make new list without the outliers
        for i in range(len(x)):
            if x[i]<limit:
                o.append(x[i])
        #avg next ten percent
        try:
            for i in range(-tenPct,0):
                a = o[i] + a
            normFactor = a/tenPct
        except IndexError:
            normFactor = 1.0
    return normFactor

def normalizeData(array, normFactor):
    newArray = [0.0]*len(array)
    for i in range(len(array)):
        # -999 indicates high background/no data
        if array[i] > -500:
            newArray[i] = array[i]/normFactor
        else:
            newArray[i] = -999.0
    return newArray

def metricAbbreviate(num):
    suffixes = {3:'k',
                6:'M',
                9:"G"}
    s = str(num)
    # replace trailing zeros with metric abbreviation
    zeroCount = len(s)-len(s.rstrip('0'))
    suffix = ''
    newString = str(s)
    for numZeros in sorted(suffixes.keys()):
        if numZeros <= zeroCount:
            suffix = suffixes[numZeros]
            newString = s[:-numZeros]
    newString = newString+suffix
    return newString
    
def writeFigures(num,reactivity,stdev,
                 rxRates,bgRates,dcRates,
                 rxErr,bgErr,dcErr,
                 rxDepth,bgDepth,dcDepth,
                 name,outputDir):
    errFlag = False
    try:
        import matplotlib as mp
        mp.rcParams["font.sans-serif"].insert(0,"Arial") # If matplotlib can find this font, then Illustrator
                                                         # should be able to open the correct font on PDF import
                                                         # rather than replacing it. Matplotlib generally falls
                                                         # back to Bitstream fonts, which Illustrator usually
                                                         # can't locate by default.
        mp.rcParams["font.family"] = "sans-serif"
        mp.rcParams["pdf.fonttype"] = 42 # use TrueType fonts when exporting pdfs (embeds most fonts)
        mp.rcParams['xtick.direction'] = 'out'
        mp.rcParams['ytick.direction'] = 'out'
        mp.rcParams['legend.fontsize'] = 14
        mp.rcParams['grid.color'] = ".8"
        mp.rcParams['grid.linestyle'] = '-'
        mp.rcParams['grid.linewidth'] = 1
        mp.use('Agg')

        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        sys.stderr.write("Warning: matplotlib python module was not found, so no figures were rendered for %i.\n"%(name))
        errFlag = True
    version = [int(v) for v in mp.__version__.split('.')]
    # Check that we have matplotlib version 1.2 or greater
    if version[0]<1 or (version[0]==1 and version[1]<2):
        sys.stderr.write("Warning: matplotlib version 1.2 or greater is recommended, but version %s was found, so no figures were rendered for %i.\n"%(mp.__version__,name))
        errFlag = True
    if errFlag==True:
        return False

    from numpy import percentile

    # Add a zeroeth nuc so axis numbering works correctly
    # There's probably a better way to do this
    num = [0]+num
    reactivity = [0]+reactivity
    stdev = [0]+stdev
    rxDepth = [0]+rxDepth
    bgDepth = [0]+bgDepth
    dcDepth = [0]+dcDepth
    rxRates = [0]+rxRates
    bgRates = [0]+bgRates
    dcRates = [0]+dcRates
    rxErr = [0]+rxErr
    bgErr = [0]+bgErr
    dcErr = [0]+dcErr

    orangeThresh = 0.4
    redThresh = 0.85

    grayVals = []
    grayNums = []
    grayErrs = []
    blackVals = []
    blackNums = []
    blackErrs = []
    orangeVals = []
    orangeNums = []
    orangeErrs = []
    redVals = []
    redNums = []
    redErrs = []
    for i in range(len(reactivity)):
        if reactivity[i] < -900:
            grayVals.append(-1)
            grayNums.append(num[i])
            grayErrs.append(0)
        elif reactivity[i] < orangeThresh:
            blackVals.append(reactivity[i])
            blackNums.append(num[i])
            blackErrs.append(stdev[i])
        elif reactivity[i] < redThresh:
            orangeVals.append(reactivity[i])
            orangeNums.append(num[i])
            orangeErrs.append(stdev[i])
        else:
            redVals.append(reactivity[i])
            redNums.append(num[i])
            redErrs.append(stdev[i])

    yMin, yMax = (-0.5, 4)
    leftInches = 0.9
    rightInches = 0.4
    spWidth = len(num)*0.032
    figWidth = max(7,spWidth+leftInches+rightInches)
    fig = plt.figure(figsize=(figWidth,8))
    leftPercent = leftInches/figWidth
    rightPercent = 1-rightInches/figWidth
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(313)
    ax3 = plt.subplot(312)
    plt.subplots_adjust(hspace=0.5, left=leftPercent,right=rightPercent,top=0.94)

    nearBlack = (0,0,1/255.0)

    if len(grayNums)>0:
        ax1.bar(grayNums,grayVals,
                 align="center",
                 width=1.05, color="0.80", edgecolor="0.80",linewidth=0.0)
    if len(blackNums)>0:
        ax1.bar(blackNums,blackVals,
                 align="center",
                 width=1.05, color="black", edgecolor="black",linewidth=0.0,
                 yerr=blackErrs,ecolor=nearBlack,capsize=1)
    if len(orangeNums)>0:
        ax1.bar(orangeNums,orangeVals,
                 align="center",
                 width=1.05, color="orange",edgecolor="orange",linewidth=0.0,
                 yerr=orangeErrs,ecolor=nearBlack,capsize=1)
    if len(redNums)>0:
        ax1.bar(redNums,redVals,
                 align="center",
                 width=1.05,color="red",edgecolor="red",linewidth=0.0,
                 yerr=redErrs,ecolor=nearBlack,capsize=1)
    
    ax1title = ax1.set_title(name, horizontalalignment="left", fontsize=16)
    x,y = ax1title.get_position()
    ax1title.set_position((0,y))
    ax1.set_ylim(yMin,yMax)
    ax1.set_xlim(1,len(num))
    #ax1.set_yticks(fontsize=9)

    #tickNums = range(num[0]+10,num[-1]+1,10)
    #tickPos = range(num[0]+9,num[-1],10)
    #ax1.set_xticks(tickPos,tickNums,fontsize=9,rotation=30)
    #ax1.set_xticks(fontsize=9)

    ax1.yaxis.grid(True)
    ax1.set_axisbelow(True)

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    for loc, spine in ax1.spines.items():
        if loc == 'bottom':
            spine.set_position(('outward', 6))  # move outward (down) 6 pts
            spine.set_smart_bounds(True)
    for loc, spine in ax1.spines.items():
        if loc == 'left':
            spine.set_position(('outward', 6))  # move outward (left) 6 pts
            spine.set_smart_bounds(True)

    # need to add labels after moving spines, otherwise they will disappear
    ax1xlabel = ax1.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax1xlabel.get_position()
    print str((x,y))
    ax1xlabel.set_position((0,y))
    ax1ylabel = ax1.set_ylabel("Shape Reactivity", horizontalalignment="left", fontsize=14)
    x,y = ax1ylabel.get_position()
    ax1ylabel.set_position((x,0))

    # add a SHAPE colorbar to the vertical axis
    # uses a little transformation magic to place correctly
    inv = ax1.transData.inverted()
    for loc, spine in ax1.spines.items():
        if loc == 'left':
            trans = spine.get_transform()
    pt = trans.transform_point([0,0])
    pt2 = inv.transform_point(pt)
    rectX = pt2[0]
    ptA = (0,0)
    ptB = (6,0)
    ptA2 = inv.transform_point(ptA)
    ptB2 = inv.transform_point(ptB)
    rectW = ptB2[0]-ptA2[0]    
    rect = Rectangle((rectX,-0.5), rectW, orangeThresh+0.5, facecolor="black", edgecolor="none")
    ax1.add_patch(rect)
    rect.set_clip_on(False)
    rect = Rectangle((rectX,orangeThresh), rectW, redThresh-orangeThresh, facecolor="orange", edgecolor="none")
    ax1.add_patch(rect)
    rect.set_clip_on(False)
    rect = Rectangle((rectX,redThresh), rectW, 4-redThresh, facecolor="red", edgecolor="none")
    ax1.add_patch(rect)
    rect.set_clip_on(False)
    
    ax1.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax1.get_yaxis().tick_left()

    ax1.tick_params(axis='y',which='minor',left='off')
    #ax1.tick_params(axis='x',which='minor')

    ax1.minorticks_on()

    yticks = ax1.get_yticks()
    #print str(yticks)
    strippedTicks = ["%d"%val if val==int(val) else "%s"%val for val in yticks]
    ax1.set_yticklabels(strippedTicks)
    
    for line in ax1.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax1.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax1.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    
    # put nuc sequence below axis
    fontProp = mp.font_manager.FontProperties(family = "monospace", style="normal",weight="bold",size="4.5")
    for i in range(len(seq)):
        nuc = seq[i]
        if nuc == "T":
            nuc = "U"
        colorDict = {"A":"#f20000",
                     "U":"#f28f00",
                     "G":"#00509d",
                     "C":"#00c200"}
        if nuc in colorDict.keys():
            col = colorDict[nuc]
        else:
            col = "black"
        ax1.annotate(nuc, xy=(i+1, -0.67),fontproperties=fontProp,color=col,annotation_clip=False, horizontalalignment="center")

    ax2.plot(num,rxDepth, linewidth = 1.5, color=rxColor)
    ax2.plot(num,bgDepth, linewidth = 1.5, color=bgColor)
    ax2.plot(num,dcDepth, linewidth = 1.5, color=dcColor)
    ax2.set_xlim(1,len(num))
    #ax2.legend(["+Reagent","Background","Denatured"], bbox_to_anchor=(1.1,1.1))
    leg = ax2.legend(["Modified","Untreated","Denatured"], loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    xmin, xmax, ymin, ymax = ax2.axis()
    ax2.set_ylim(0,ymax)
    #ax2.set_yscale('symlog')# useful but currently disabled because of a matplotlib/pyparsing bug
    ax2xlabel = ax2.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax2xlabel.get_position()
    ax2xlabel.set_position((0,y))
    ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=14)
    x,y = ax2ylabel.get_position()
    ax2ylabel.set_position((x,0))
    #ax2title = ax2.set_title(name, horizontalalignment="left")
    #x,y = ax2title.get_position()
    #ax2title.set_position((0,y))

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax2.get_yaxis().tick_left()

    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='minor',left='off')
    #ax2.tick_params(axis='x',which='minor')

    #xlabels = ["%.2f"%v for v in xticks]
    #ax3.set_xticks(xticks)
    #ax3.set_xticklabels(xlabels,rotation = -45, horizontalalignment='left')

    yticks = [int(y) for y in ax2.get_yticks()]
    formattedTicks = []
    for val in yticks:
        formattedTicks.append(metricAbbreviate(val))
    ax2.set_yticklabels(formattedTicks)

    for line in ax2.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax2.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax2.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax2.yaxis.grid(True)
    ax2.set_axisbelow(True)

    # choose a decent range for axis, excluding high-background positions
    goodIndices = []
    for i in xrange(len(bgRates)):
        if bgRates[i]<=conf.maxBackground:
            goodIndices.append(i)
    tempRates = [rxRates[i] for i in goodIndices]
    nearTopRate = percentile(tempRates,98.0)
    maxes = [0.32,0.16,0.08,0.04,0.02,0.01]
    yMax = maxes[0]
    for i in xrange(len(maxes)):
        if nearTopRate<maxes[i]:
            yMax = maxes[i]

    rxUpper = [rxRates[i]+rxErr[i] for i in xrange(len(rxRates))]
    rxLower = [rxRates[i]-rxErr[i] for i in xrange(len(rxRates))]
    bgUpper = [bgRates[i]+bgErr[i] for i in xrange(len(bgRates))]
    bgLower = [bgRates[i]-bgErr[i] for i in xrange(len(bgRates))]
    dcUpper = [dcRates[i]+dcErr[i] for i in xrange(len(dcRates))]
    dcLower = [dcRates[i]-dcErr[i] for i in xrange(len(dcRates))]

    ax3xlabel = ax3.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax3xlabel.get_position()
    ax3xlabel.set_position((0,y))
    ax3ylabel = ax3.set_ylabel("Mutation rate (%)", horizontalalignment="left", fontsize=14)
    x,y = ax3ylabel.get_position()
    ax3ylabel.set_position((x,0))
    
    ax3.plot(rxRates, zorder=3,color=rxColor,linewidth=1.5)
    ax3.plot(bgRates, zorder=2, color=bgColor, linewidth=1.5)
    ax3.fill_between(num,rxLower,rxUpper,edgecolor="none",alpha=0.5, facecolor=rxColor)
    ax3.fill_between(num,bgLower,bgUpper,edgecolor="none",alpha=0.5, facecolor=bgColor)
    ax3.legend(["Modified","Untreated"], loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    ax3.set_xlim((1,len(rxRates)))
    ax3.set_ylim((0,yMax))

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.minorticks_on()
    ax3.tick_params(axis='y',which='minor',left='off')

    ticks = [x*100 for x in ax3.get_yticks()]
    ax3.set_yticklabels(["%d"%val if val==int(val) else "%s"%val for val in ticks])

    for line in ax3.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax3.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax3.yaxis.grid(True)
    ax3.set_axisbelow(True)

    # TODO: add a tick for the first nuc - can't seem to add one without screwing
    # up all the other ticks

    outPath = os.path.join(outputDir,"reactivity_profiles")
    outPath = os.path.join(outPath,name+"_depth_and_reactivity.pdf")
    plt.savefig(outPath)
    return True

def commify(s):
    o = ""
    for i in xrange(-1,-len(s)-1,-1):
        if -i%3==0:
            o = ','+s[i]+o
        else:
            o = s[i]+o
    o = o.lstrip(',')
    return o 

def drawMedian(ax,vals,col, intOnly=False):
    xmin, xmax, ymin, ymax = ax.axis()
    med = numpy.median(vals)
    ax.axvline(x=med, color=col,zorder=0,alpha=0.7)
    txt = str(med)
    if intOnly==True:
        txt = metricAbbreviate(int(med))
    ax.annotate(txt,xy=(med,(ymax-ymin)*0.6), color=col, rotation=-45, fontweight="bold")

def drawPercentile(ax,vals,percent, col, x1=0, x2=0, y=0, percentage=False, name=""):
    from numpy import percentile
    xmin, xmax, ymin, ymax = ax.axis()
    #med = calcQuartile(vals,percentile)
    med = percentile(vals, percent)
    if percentage==True:
        med = med*100
    if med==int(med):
        txt="%d"%med
    else:
        txt="%0.2f"%med
        if percentage==True:
            txt=txt.rstrip('0')
    if percentage==False:
        txt=commify("%i"%med)
    if percentage==True:
        txt = txt+"%"
    ax.text(x1, y, name, transform=ax.transAxes, fontsize=10, horizontalalignment="left")
    ax.text(x2, y, txt, transform=ax.transAxes, fontsize=10, horizontalalignment="right")


def writeHistograms(num,normShape,normStderr,
                    rxRate,bgRate,dcRate,
                    rxDepth,bgDepth,dcDepth,
                    name,outputDir):
    from numpy import percentile
    from math import ceil
    import matplotlib as mp
    mp.rcParams["font.sans-serif"].insert(0,"Arial")
    mp.rcParams["font.family"] = "sans-serif"
    mp.rcParams["pdf.fonttype"] = 42 # use TrueType fonts when exporting pdfs (embeds most fonts)
    mp.rcParams['xtick.direction'] = 'out'
    mp.rcParams['ytick.direction'] = 'out'
    mp.rcParams['legend.fontsize'] = 14
    mp.rcParams['grid.color'] = ".8"
    mp.rcParams['grid.linestyle'] = '-'
    mp.rcParams['grid.linewidth'] = 1
    mp.use('Agg')

    import matplotlib.pyplot as plot

    tLabelSize = 10

    filterIndices = [i for i in xrange(len(normShape)) if normShape[i] > -900]
    filteredBgRate = [bgRate[i] for i in filterIndices]
    filteredRxRate = [rxRate[i] for i in filterIndices]
    filteredDcRate = [dcRate[i] for i in filterIndices]

    filteredShape = [normShape[i] for i in filterIndices]

    filteredStderr = [normStderr[i] for i in filterIndices]

    fig, (ax1,ax2,ax3) = plot.subplots(nrows=1, ncols=3)
    fig.set_size_inches(10,6)
    fig.subplots_adjust(bottom=0.6, top=0.85, wspace=0.5, left=0.08, right=0.95)

    plot.suptitle(name,fontsize=18,horizontalalignment="left", x=0.02)
    
    #max90percentile = max([calcQuartile(filteredBgRate,0.90),
    #                       calcQuartile(filteredRxRate,0.90),
    #                       calcQuartile(filteredDcRate,0.90)])
    max90percentile = max([percentile(filteredBgRate,90.0),
                           percentile(filteredRxRate,90.0),
                           percentile(filteredDcRate,90.0)])

    numBins = 30
    # use fewer bins for short RNAs
    if len(filteredShape)<500:
        numBins = 10
    intBins = range(0,numBins+1,1)
    
    # adjust rate axis for higher mutation rates (e.g. from DMS or other highly mutagenic reagents)
    binMax = 0.03
    if max90percentile < 0.01:
        binMax = 0.008
    binWidth = binMax/(numBins)
    floatBins = [b*binWidth for b in intBins]
    
    ax1.set_ylabel("Normalized\nnucleotide count",fontsize=13)
    ax1.set_xlabel("Mutation rate (%)", fontsize=13)
    #title = ax1.set_title("Mutation rates",fontsize=15)
    ax1.set_title('Mutation rates', x=0.5,y=1.08) 
    #y = title.get_y()
    #title.set_y(y

    ax1.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax1.get_yaxis().tick_left()

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    for line in ax1.get_yticklines() + ax1.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)
    
    ax1.tick_params(axis='both', labelsize=tLabelSize)
    
    rxpdf, bins, rxRatePatches = ax1.hist(filteredRxRate, floatBins, histtype='step', lw=1.5,color=rxColor,zorder=4)
    bgpdf,bins, bgRatePatches = ax1.hist(filteredBgRate, floatBins, histtype='step', lw=1.5, color=bgColor,zorder=3)
    dcpdf,bins, dcRatePatches = ax1.hist(filteredDcRate, floatBins, histtype='step', lw=1.5, color=dcColor,zorder=2)    

    maxRxRateY = max([x[1] for x in rxRatePatches[0].get_xy()])
    maxBgRateY = max([x[1] for x in bgRatePatches[0].get_xy()])
    maxDcRateY = max([x[1] for x in dcRatePatches[0].get_xy()])

    # rescale steplots to max out at 1
    rxRatePatches[0].set_xy([[c[0],c[1]/maxRxRateY] for c in rxRatePatches[0].get_xy()])
    bgRatePatches[0].set_xy([[c[0],c[1]/maxBgRateY] for c in bgRatePatches[0].get_xy()])
    dcRatePatches[0].set_xy([[c[0],c[1]/maxDcRateY] for c in dcRatePatches[0].get_xy()])
    
    ymax = 1.0
    xmax = binMax
    
    leg = ax1.legend(["Modified","Untreated","Denatured"], framealpha=0.75, fontsize=11)
    for l in leg.get_lines():
        l.set_linewidth(1.5)
        l.set_alpha(1.0)
    
    ax1.set_xlim([0.0,xmax])
    ax1.set_ylim([0,ymax])

    ticks = [x*100 for x in ax1.get_xticks()]
    ax1.set_xticklabels(["%d"%val if val==int(val) else "%s"%val for val in ticks])

    drawPercentile(ax1,filteredRxRate,95.0, "black", x1=0, x2=1, y=-1, percentage=True, name="Modified sample:\n 95th percentile rate:")
    drawPercentile(ax1,filteredRxRate,5.0, "black", x1=0,x2=1, y=-1.1, percentage=True, name=" Median rate:")

    drawPercentile(ax1,filteredBgRate,95.0, "black", x1=0,x2=1, y=-1.4, percentage=True, name="Untreated sample:\n 95th percentile rate:")
    drawPercentile(ax1,filteredBgRate,5.0, "black", x1=0,x2=1, y=-1.5, percentage=True, name=" Median rate:")

    drawPercentile(ax1,filteredDcRate,95.0, "black", x1=0,x2=1, y=-1.8, percentage=True, name="Denatured sample:\n 95th percentile rate:")
    drawPercentile(ax1,filteredDcRate,5.0, "black", x1=0,x2=1, y=-1.9, percentage=True, name=" Median rate:")

    ax2.set_title("Read depths", x=0.5,y=1.08)
    ax2.set_ylabel("Normalized\nnucleotide count",fontsize=13)
    ax2.set_xlabel("Read depth", fontsize=13)

    max90percentile = max([percentile(bgDepth,90.0),
                           percentile(rxDepth,90.0),
                           percentile(dcDepth,90.0)])
    
    xmax = max90percentile*1.1
    
    rxpdf, bins,rxDepthPatches = ax2.hist(rxDepth, bins=numBins, range=(0,xmax) ,histtype='step', lw=1.5,color=rxColor,zorder=4)
    bgpdf,bins,bgDepthPatches =  ax2.hist(bgDepth, bins=numBins, range=(0,xmax) ,histtype='step', lw=1.5,color=bgColor,zorder=3)
    dcpdf, bins,dcDepthPatches = ax2.hist(dcDepth, bins=numBins, range=(0,xmax), histtype='step',lw=1.5,color=dcColor,zorder=2)

    maxRxDepthY = max([x[1] for x in rxDepthPatches[0].get_xy()])
    maxBgDepthY = max([x[1] for x in bgDepthPatches[0].get_xy()])
    maxDcDepthY = max([x[1] for x in dcDepthPatches[0].get_xy()])

    # rescale steplots to max out at 1
    rxDepthPatches[0].set_xy([[c[0],c[1]/maxRxDepthY] for c in rxDepthPatches[0].get_xy()])
    bgDepthPatches[0].set_xy([[c[0],c[1]/maxBgDepthY] for c in bgDepthPatches[0].get_xy()])
    dcDepthPatches[0].set_xy([[c[0],c[1]/maxDcDepthY] for c in dcDepthPatches[0].get_xy()])
    
    ymax = 1.0
                           
    #leg = ax2.legend(["Modified","Untreated","Denatured"], framealpha=0.75, fontsize=11)
    #for l in leg.get_lines():
    #    l.set_linewidth(1.5)
    #    l.set_alpha(1.0)

    ax2.set_xlim([0.0,xmax])
    ax2.set_ylim([0.0,ymax])
    
    drawPercentile(ax2,rxDepth,50.0, "black", x1=0,x2=1, y=-1, name="Modified sample:\n Median depth:")
    drawPercentile(ax2,rxDepth,5.0, "black", x1=0,x2=1, y=-1.1, name=" 5th percentile depth:")

    drawPercentile(ax2,bgDepth,50.0, "black", x1=0,x2=1, y=-1.4, name="Untreated sample:\n Median depth:")
    drawPercentile(ax2,bgDepth,5.0, "black", x1=0,x2=1, y=-1.5, name=" 5th percentile depth:")

    drawPercentile(ax2,dcDepth,50.0, "black", x1=0,x2=1, y=-1.8, name="Denatured sample:\n Median depth:")
    drawPercentile(ax2,dcDepth,5.0, "black", x1=0,x2=1, y=-1.9, name=" 5th percentile depth:")

    ax2.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax2.get_yaxis().tick_left()

    ax2.tick_params(axis='both', labelsize=tLabelSize) 

    xticks = [int(x) for x in ax2.get_xticks()]
    formattedTicks = []
    for val in xticks:
        formattedTicks.append(metricAbbreviate(val))
    ax2.set_xticklabels(formattedTicks, rotation=90)

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    for line in ax2.get_yticklines() + ax2.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)

    ax3.set_xlabel("Normalized reactivity", fontsize=13)
    ax3.set_ylabel("Normalized\nnucleotide count", fontsize=13)
    ax3.set_title("Reactivity distribution", x=0.5,y=1.08)
    #xticks = [x/10.0 for x in range(-5, 30, 5)] + [0.4, 0.85]
    #xticks = [-1.0, -0.5, 0.0, 0.4, 0.85, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    xticks = [-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    xlabels = ["%d"%val if val==int(val) else "%s"%val for val in xticks]

    ax3.tick_params(axis='both', labelsize=tLabelSize) 

    ax3.set_xticks(xticks)
    ax3.set_xticklabels(xlabels)

    pdf, bins,shapePatches = ax3.hist(filteredShape, bins=numBins, range=(-1,4), histtype='step', lw=1.5, color="black", zorder=4)
    pdf, bins, errPatches = ax3.hist(filteredStderr, bins=numBins, range=(-1,4), histtype='step', lw=1.5, color='blue', zorder=3)

    maxShapeY = max([x[1] for x in shapePatches[0].get_xy()])
    maxErrY = max([x[1] for x in errPatches[0].get_xy()])

    # rescale steplots to max out at 1
    shapePatches[0].set_xy([[c[0],c[1]/maxShapeY] for c in shapePatches[0].get_xy()])
    errPatches[0].set_xy([[c[0],c[1]/maxErrY] for c in errPatches[0].get_xy()])

    leg = ax3.legend(["Reactivities","Stderrs"], framealpha=0.75, fontsize=11)

    ax3.axvline(x=0,ls='-',lw=1,color="0.5",zorder=0)
    #ax3.axvline(x=0.4,ls='-',lw=1,color="0.5",zorder=0)
    #ax3.axvline(x=0.85,ls='-',lw=1,color="0.5",zorder=0)

    #ymin, ymax = ax3.get_ylim()
    #ymax = max(shapePdf)
    ymax = 1
    #ax3.set_xlim([-0.5,3])
    ax3.set_ylim([0,ymax])

    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    for line in ax3.get_yticklines() + ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(1)
    
    outPath = os.path.join(outputDir,"reactivity_profiles")
    outPath = os.path.join(outPath,name+"_histograms.pdf")
    plot.savefig(outPath)
    
    return
    
def percentile(N, P):
## modified from http://code.activestate.com/recipes/511478/
# Find the percentile of a list of values
# N - A list of values.  N must be sorted.
# P - A float value from 0.0 to 1.0
    n = int(round(P * len(N) + 0.5))
    val = N[n-1]
    return val

def filterProfile(profile, stderr,bgSignal,rxDepth,bgDepth,dcDepth,
                  maxBackground=0.05, minDepth=0):
# filter reactivity profile -
# exclude nucs with high background mutation rates
# exclude nucs with low read coverage
# exclude nucs with high standard errors (optional)

    filteredProfile = [0.0]*len(profile)
    filteredStderr = [0.0]*len(profile)

    for i in range(len(profile)):
        goodNuc = True
        if bgSignal[i] > maxBackground:
            goodNuc = False
        if bgDepth[i] < minDepth or rxDepth[i] < minDepth or dcDepth[i] < minDepth:
            goodNuc = False
        if conf.filterStderr == True:
            if (stderr[i]>(abs(profile[i])*0.5+0.4)):
                goodNuc = False

        if profile[i] > -500 and goodNuc == True:
            filteredProfile[i] = profile[i]
            filteredStderr[i] = stderr[i]
        else:
            filteredProfile[i] = -999.0
            filteredStderr[i] = 0.0
    return filteredProfile, filteredStderr

def safeDivide(num, den):
    # divide each element in a list by the corresponding element in another list
    # if this fails, report -999
    outList = [-999.0]*len(num)
    for i in range(len(num)):
        if num[i] > -500 and den[i] > -500:
            try:
                outList[i] = num[i]/float(den[i])
            except ZeroDivisionError:
                pass
    return outList

try:
    import conf    
    
    csvDir = sys.argv[1]
    profileName = sys.argv[2]
    targetName = sys.argv[3]
    rxNames = sys.argv[findArg("--plus",sys.argv)+1:findArg("--minus",sys.argv)]
    bgNames = sys.argv[findArg("--minus",sys.argv)+1:findArg("--denat",sys.argv)]
    dcNames = sys.argv[findArg("--denat",sys.argv)+1:]
    
    #print "profileName: %s \n"%profileName
    #print "targetName: %s \n"%targetName
    #print "rxNames:"
    #for name in rxNames:
    #    print name+"\n"
    #print "minusNames:"
    #for name in bgNames:
    #    print name+"\n"
    #print "denatNames:"
    #for name in dcNames:
    #    print name+"\n"
    

    csvFileList = os.listdir(csvDir)
    csvFileList = [filename for filename in csvFileList if os.path.splitext(filename)[1]==".csv"]
    csvFilePathList = [os.path.join(csvDir,filename) for filename in csvFileList]
    # {sample name:file path}
    #sampleDict = {}
    #for filename in csvFileList:
    #    sampleDict[os.path.splitext(filename)[0]] = os.path.join(csvDir,filename)

    # make combined plus reagent mutation freq list    
    #print "sampleDict: "+str(sampleDict)
    rxRates,rxDepth,seq = combineData(csvFileList, csvFilePathList, rxNames)

    # make combined minus reagent mutation freq list
    bgRates, bgDepth, unused = combineData(csvFileList, csvFilePathList, bgNames)

    # make combined denatured control mutation freq list
    dcRates, dcDepth, unused = combineData(csvFileList, csvFilePathList, dcNames)
        
    # create reactivity profile
    diff = [rxRates[i]-bgRates[i] for i in range(len(rxRates))]     
    profile = safeDivide(diff,dcRates)
     
    # create error estimates for final profile
    stderr = [0.0]*len(rxRates)
    rxStderr = [0.0]*len(rxRates)
    bgStderr = [0.0]*len(rxRates)
    dcStderr = [0.0]*len(rxRates)
    for i in range(len(rxRates)):
        # these are estimates of the standard error of the mean for each condition
        # Since the data are Poisson distributed, sample mean (i.e. mutation rate) = variance
        try:
            rxStderr[i] = math.sqrt(rxRates[i])/math.sqrt(float(rxDepth[i]))
            bgStderr[i] = math.sqrt(bgRates[i])/math.sqrt(float(bgDepth[i]))
            dcStderr[i] = math.sqrt(dcRates[i])/math.sqrt(float(dcDepth[i]))
        except ValueError:
            stderr[i] = 0.0
            continue
        
        try:
            stderr[i] = math.sqrt( (rxStderr[i]/dcRates[i])**2 + (bgStderr[i]/dcRates[i])**2 + (dcStderr[i]*(rxRates[i]-bgRates[i])/(dcRates[i]**2))**2 )
        except ZeroDivisionError:
            stderr[i] = 0.0
            continue

    # handle excluded values
    for i in range(len(rxRates)):
        if rxRates[i] < -900 or bgRates[i] < -900 or dcRates[i] < -900:
            profile[i] = -999.0
            stderr[i] = 0.0

    #print "\t".join([str(n) for n in plusFreqs])
    #print "\t".join([str(n) for n in minusFreqs])
    #print "\t".join([str(n) for n in denatFreqs])

    # filter nucs with high background or low denatured control
    filteredProfile, filteredStderr = filterProfile(profile, stderr, bgRates, rxDepth, bgDepth, dcDepth,
                                    maxBackground=conf.maxBackground, minDepth=conf.minDepth)

    #print ",".join([str(n) for n in profile])
    #print ",".join([str(n) for n in filteredProfile])        

    # normalize profile to standard scale
    # exclude low-coverage regions from calculation of norm factor
    tempProfile = [filteredProfile[i] for i in xrange(len(filteredProfile)) if min([bgDepth[i],rxDepth[i],dcDepth[i]])>conf.minNormDepth]
    # warn if fewer than 60 nucleotides remain for normalization calculation
    if len(tempProfile) < 60:
        sys.stderr.write("Warning: only %i nucleotides are above the minNormDepth threshold of %i in reactivity profile %s. Normalization of this profile may be unreliable.\n"%(len(tempProfile),conf.minNormDepth,profileName))
    normFactor = findBoxplotFactor(tempProfile)
    if conf.normProfile == False:
        normFactor = 1.0
    normProfile = normalizeData(filteredProfile, normFactor) 

    # adjust estimated stdevs by norm factor
    normStderr = [0.0]*len(stderr)
    for i in range(len(stderr)):
        if normProfile[i] < -900:
            normStderr[i] = 0.0
        else:
            normStderr[i] = stderr[i]/normFactor

    sortedStderr = sorted([normStderr[i] for i in range(len(normStderr)) if normProfile[i]>-900])

    normedMedianStderr = numpy.median(sortedStderr)
    normedStdevStderr = numpy.std(sortedStderr)

    upperQuartile = percentile(sortedStderr, 0.75)
    lowerQuartile = percentile(sortedStderr, 0.25)

    IQRStderr = upperQuartile-lowerQuartile

    #print ",".join([str(n) for n in profile])
    #print ",".join([str(n) for n in filteredProfile])
    #print ",".join([str(n) for n in normProfile])

    # write .SHAPE file
    fileOutPath = os.path.join(conf.outputDir,"reactivity_profiles/"+profileName+".shape")
    fileOut = open(fileOutPath, "w")
    nuc = 1
    for value in normProfile:
        fileOut.write("%i\t%f\n"%(nuc,value))
        nuc += 1
    fileOut.close()

    # write .tab file
    fileOutPath = os.path.join(conf.outputDir,"reactivity_profiles/"+profileName+".tab")
    fileOut = open(fileOutPath, "w")
    fileOut.write("Nucleotide\tSequence\tbg depth\tbg rate\t")
    fileOut.write("rx depth\trx rate\tdc depth\tdc rate\trx-bg\t")
    fileOut.write("(rx-bg)/dc\tstderr\tFiltered Reactivity\tFiltered StdErr\tNormalized Reactivity\tNormalized StdErr\tMedian Normed Stderr\tIQR of Normed Stderr\tStdev of Normed Stderr\t\n")
    nuc = 1
    for i in range(len(normProfile)):
        #fileOut.write("%i,%s,%f,%f\n"%(nuc,seq[i],normProfile[i],normStdev[i]))
        fileOut.write("%i\t%s\t%i\t%f\t%i\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t\n"%(
            nuc,seq[i],
            bgDepth[i],bgRates[i],
            rxDepth[i],rxRates[i],
            dcDepth[i],dcRates[i],
            diff[i],
            profile[i],stderr[i],
            filteredProfile[i],filteredStderr[i],
            normProfile[i],normStderr[i],
            str(normedMedianStderr) if i==0 else "",
            str(IQRStderr) if i==0 else "",
            str(normedStdevStderr) if i==0 else ""))
        nuc += 1
    fileOut.close()

    # write .map file
    fileOutPath = os.path.join(conf.outputDir, "reactivity_profiles/"+profileName+".map")
    fileOut = open(fileOutPath, "w")
    nuc = 1
    for i in range(len(normProfile)):
        fileOut.write("%i\t%f\t%f\t%s\n"%(nuc,normProfile[i],normStderr[i],seq[i]))
        nuc += 1
    fileOut.close()

    # write reactivity profile and coverage depth figures    
    num = range(1,len(normProfile)+1)
    success= writeFigures(num,normProfile,normStderr,
                          rxRates,bgRates,dcRates,
                          rxStderr,bgStderr,dcStderr,
                          rxDepth,bgDepth,dcDepth,
                          profileName,conf.outputDir)
    if success==True:
        writeHistograms(num,normProfile,normStderr,
                        rxRates,bgRates,dcRates,
                        rxDepth,bgDepth,dcDepth,
                        profileName,conf.outputDir)

except Exception as e:
    sys.stderr.write(traceback.format_exc())
    exit(1)
