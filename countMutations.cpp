// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// Counts mutations and sequencing depths.
// Copyright Steven Busan 2014

/*---------------------------------------------------------------------
GPL statement:

This file is part of Shapemapper.

ShapeMapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ShapeMapper is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ShapeMapper.  If not, see <http://www.gnu.org/licenses/>.
*///-------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iterator>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include "string_funcs.h"
#include "ref_seq.h"

using namespace std;


struct ParsedReadFields
{
    string clusterName;
    int nucStart;
    int nucEnd;
    string seqCodesPruned;
    string seqCodesUnpruned;
    vector<int> quals;
};

class ParsedReadInputFile
{
    // encapsulates simple file operations on a parsed read file
public:
    ParsedReadInputFile(const string& filePath);
    ifstream fileIn;
    int readLine(string& line);
    int splitFields(const string& line, ParsedReadFields& fields);
};

ParsedReadInputFile::ParsedReadInputFile(const string& filePath)
{
    fileIn.open(filePath.c_str());
    if (fileIn.is_open())
    {
        return;
    }else
    { 
	std::stringstream message;
        message << "Unable to open file " << filePath << std::endl << std::flush;
        throw std::runtime_error(message.str().c_str()); 
    }
}

int ParsedReadInputFile::readLine(string& line)
{
    // get a line from the file, store in string passed as argument
    // if there are no more lines in the file, set the end-of-file return flag to 1
    int eofFlag = 0;
    getline( fileIn, line);
    if(!fileIn.good())
    { 
        eofFlag = 1;
        return eofFlag;
    }else return eofFlag;
}

int ParsedReadInputFile::splitFields(const string& line, ParsedReadFields& fields)
{
    vector<string> splitLine;
    splitLineTab(line, splitLine);
    //cout << "line split by tabs\n" << flush;
    //cout << "splitLine.size() = " << splitLine.size() << "\n" << flush;
    if (splitLine.size() < 6){ 
        //cout << line << "\n" <<flush;
    }
    fields.clusterName = splitLine[0];
    fields.nucStart = atoi(splitLine[1].c_str());
    fields.nucEnd = atoi(splitLine[2].c_str());
    fields.seqCodesPruned = splitLine[3];
    fields.seqCodesUnpruned = splitLine[4];
    fields.quals.resize(splitLine[5].size());
    int i;
    for (i=0; i<splitLine[5].size(); i++)
    {
        fields.quals[i] = (int)(splitLine[5][i])-33;
    }
    return 0;
}

float safeDivide(const int& count, const int& divisor)
{
    if (divisor==0) return -999.0;
    else return (float)count/(float)divisor;
}

class MutationCounters
{
public:
    MutationCounters(int length);
    vector<int> depth;
    map<char,map<char, vector<int> > > mismatchCount;
    map<char,vector<int> > deletionCount;
    vector<int> mutationCount;
    vector<float> mutationRate;
    int totalMismatches;
    int totalDeletions;
    int phredRejectedMutations;
    int calcMutationRates();
private:
    int charsToVectorIndex(char fromNuc, char toNuc);
};

MutationCounters::MutationCounters(int length)
{
    depth.resize(length);
    mutationCount.resize(length);
    mutationRate.resize(length);
    totalMismatches = 0;
    totalDeletions = 0;
    phredRejectedMutations = 0;
    char nucs [4] = {'A','T','G','C'};
    char fromNuc;
    char toNuc;
    int i;
    int j;
    for(i=0; i<4; i++)
    {
        fromNuc = nucs[i];
        for(j=0; j<4; j++)
        {
            deletionCount[fromNuc].resize(length);
            toNuc = nucs[j];
            //if (fromNuc != toNuc)
            //{
        
            mismatchCount[fromNuc][toNuc].resize(length);
	    
            //}
        }
    }
}

int MutationCounters::calcMutationRates()
{
    int i;
    for(i=0; i<depth.size(); i++)
    {
        mutationRate[i] = safeDivide(mutationCount[i],depth[i]);
    }
    return 0;
}

bool inBounds(const int& i, const string& seq)
{
    if (i>=0 and i<seq.length()) return 1;
    else return 0;
}

int countMutationsInRead(const int startNuc,
                         const string& seqCodes,
                         const string& targetSeq,
                         const vector<int>& quals,
                         const int& minPhred, 
                         MutationCounters& mutCount)
{
    int i;
    int refIndex;
    int q;
    char c,d;
    //cout << "seq length = " << targetSeq.length() << "\n" << flush;
    for (i=0; i<seqCodes.length()-1; i++)
    {
        refIndex = i+startNuc-1; // 0-based
        c = seqCodes[i];
        q = quals[i];
        d = seqCodes[i+1];
        if (c!='s' and c!='~' and q>=minPhred)
        {
            if (inBounds(refIndex,targetSeq)){ mutCount.depth[refIndex]++; }
        }
        if ((c=='A' or c=='T' or c=='G' or c=='C') and 
            (d!='A' and d!='T' and d!='G' and d!='C' and d!='-'))
        {
            //if (c==targetSeq[refIndex]) cout << "mismatch not a mismatch: nuc " << c << ", refIndex " << refIndex<<"\n"<<flush;
	    if (q>=minPhred){
                if (inBounds(refIndex,targetSeq)){ 
                    mutCount.mismatchCount[targetSeq[refIndex]][c][refIndex]++;
                    mutCount.mutationCount[refIndex]++;
                    mutCount.totalMismatches++;
		}
            }else{ mutCount.phredRejectedMutations++; }
        }
        else if (c=='-' and 
                 d!='A' and d!='T' and d!='G' and d!='C' and d!='-')
        {
	    if (q>=minPhred){
                if (inBounds(refIndex,targetSeq)){
                    //cout << "deleted nuc = "<<targetSeq[refIndex]<< ", refIndex = "<<refIndex<<"\n"<<flush;
                    mutCount.deletionCount[targetSeq[refIndex]][refIndex]++;
                    mutCount.mutationCount[refIndex]++;
                    mutCount.totalDeletions++;
	        }
            }else{ mutCount.phredRejectedMutations++; }
        }
    }
    mutCount.calcMutationRates();
    return 0;
}
                              

int writeCSV(string& filePath,
             string& sampleName,
             string& targetName, 
             RefSeq& refSeq, 
             MutationCounters& unprunedCounts,
             MutationCounters& prunedCounts)
{
    char nucs [4] = {'A','T','G','C'};
    char fromNuc;
    char toNuc;
    int i;
    int j;
    int n;
    
    int seqLen = refSeq.seqs[targetName].length();

    ofstream f(filePath.c_str());
    if(!f.is_open())
    {
	std::stringstream message;
        message << "Unable to open file " << filePath << std::endl << std::flush;
        throw std::runtime_error(message.str().c_str()); 
    }

    // write nuc numbers
    f << sampleName << ',' << targetName << ",number,";
    for(n=0; n<seqLen; n++)
    {
        f << (n+1) << ',';
    }
    f << '\n' << flush;
    
    // write sequence
    f << sampleName << ',' << targetName << ",sequence," ;
    for(n=0; n<seqLen; n++)
    {
        f << refSeq.seqs[targetName][n] << ',';
    }
    f << '\n' << flush;

    // write deletion counts
    for(i=0; i<4; i++)
    {
        fromNuc = nucs[i];
        f << sampleName << ',' << targetName << ',' << "del " << fromNuc << ',';
        for(n=0; n<seqLen; n++)
        {
            f << prunedCounts.deletionCount[fromNuc][n] << ',';
        }
        f << '\n' << flush;
    }

    // write mismatch counts
    for(i=0; i<4; i++)
    {
        fromNuc = nucs[i];
        for(j=0; j<4; j++)
        {   
            toNuc = nucs[j];
            if (fromNuc != toNuc)
            {
                f << sampleName << "," << targetName << "," << fromNuc <<"->"<< toNuc << ',';
                for(n=0; n<seqLen; n++)
                {
                    f << prunedCounts.mismatchCount[fromNuc][toNuc][n] << ',';
		}
	        f << '\n' << flush;
            }
        }
    }

    //f.setf( ios::fixed, ios::floatfield );
    f.precision(14); // set the number of digits after the decimal to output    

    // write mutation rates
    f << sampleName << ',' << targetName << ",mutation rate,";
    for(n=0; n<seqLen; n++)
    {
        f <<prunedCounts.mutationRate[n] << ',';
    }
    f << '\n' << flush;

    // write coverage depths
    f << sampleName << ',' << targetName << ",depth,";
    for(n=0; n<seqLen; n++)
    {
        f << prunedCounts.depth[n] << ',';
    }
    f << '\n' << flush;

    // write total mismatches
    f << sampleName << ',' << targetName << ",total mismatches," << prunedCounts.totalMismatches << ",\n" << flush;

    // write total deletions
    f << sampleName << ',' << targetName << ",total deletions," << unprunedCounts.totalDeletions << ",\n" << flush;

    // write unambiguous deletions
    f << sampleName << ',' << targetName << ",unambiguous deletions," << prunedCounts.totalDeletions << ",\n" << flush;

    // write count of mutations failing to pass phred score cutoff
    f << sampleName << ',' << targetName << ",mutations failing phred cutoff," << prunedCounts.phredRejectedMutations << ",\n" << flush;

    f<<"\n\n"<<flush;
    
    return 0;
}

static void printUsage()
{
    cerr<<"Usage example:\n    countMutations -file_in sample1_target1.txt"<< 
          " -ref_seqs test.fa -sample_name sample1 -target_name target1"<<
          " -file_out sample1_target1_counted.csv \n"<<flush;
}

int main(int argc, char* argv[])
{
    // test phred score conversion
    /*char c = '#';
    int scoreInt = (int)c-33;
    unsigned scoreUnsigned = (unsigned)c-33;
    cout << "int score: "<<scoreInt<<"\n";
    cout << "unsigned score: "<<scoreUnsigned<<"\n"<<flush;
    exit(0);*/

    // parse commandline arguments ---------------------------
    string file_in_path = "";
    string ref_seqs_path = "";
    string sample_name = "";
    string target_name = "";
    string file_out_path = "";
    int min_phred = 0;
    
    char* endptr;
    int base = 10;    
    int optionIndex = 0;
    static struct option longOptions[] = {
	{"file_in",      required_argument, NULL, 'f' },
	{"ref_seqs",     required_argument, NULL, 's' },
	{"sample_name", required_argument, NULL, 'n' },
	{"target_name",    required_argument, NULL, 't' },
	{"file_out", required_argument, NULL, 'o' },
        {"min_phred", required_argument, NULL, 'p'},
	{NULL,    no_argument,            NULL,  0 }
    };
    static const char *optString = ":f:s:n:t:o";

    static int flags[6] = { 0 };

    while (getopt_long_only(argc, argv, optString, longOptions, &optionIndex) != -1)
    {
        switch(optionIndex)
        {
	    case 0:
            flags[0] = 1;
            file_in_path = string(optarg);
            break;
        case 1:
            flags[1] = 1;
            ref_seqs_path = string(optarg);
            break;
        case 2:
            flags[2] = 1;
            sample_name = string(optarg);
            break;
        case 3:
            flags[3] = 1;
            target_name = string(optarg);
            break;
        case 4:
            flags[4] = 1;
            file_out_path = string(optarg);
            break;
        case 5:
	    flags[5] = 1;
            errno = 0;
            min_phred = strtol(optarg, &endptr, base);
            if(errno > 0){ cerr<<"Unable to parse primer_length argument.\n"<<flush; exit(1); }
            break;
        case 'h':
        case '?':
        default:
            printUsage();
	        exit(0);
            break;
        }
    }
    bool foundAllFlags = true;
    size_t i;
    for (i=0; i<6; i++)
    {
        if (!flags[i]){ foundAllFlags = false; }
    }
    if (!foundAllFlags){
        cerr << "All arguments are required.\n"<<flush;
        printUsage();
        exit(1);
    }

    //cout << "minPhred: "<<min_phred<<"\n"<<flush;
    //exit(0);
    // --------------------------------------------------------
    //cout << "started.\n" << flush;
    
    RefSeq refSeq(ref_seqs_path);
    //refSeq.printSequences();
    //refSeq.seqs[targetName]

    //cout << "ref seqs loaded.\n" << flush;
    
    ParsedReadInputFile parsedReadFile(file_in_path);
    //cout << "starting read input\n" << flush;

    MutationCounters unprunedCounts(refSeq.seqs[target_name].length());
    MutationCounters prunedCounts(refSeq.seqs[target_name].length());

    //cout << "mutation counting data structures initialized\n" << flush;

    int lineCount = 0;
    string line;
    ParsedReadFields fields;
    while( parsedReadFile.readLine(line) == 0){
        //cout << "got line\n" << flush;
        lineCount++;
        //if (isOnlyWhitespace(line)){ cout<<"blank line at linecount "<<lineCount<<"\n"<<flush; continue;}
        parsedReadFile.splitFields(line, fields);
        //cout << "split fields\n" << flush;
        // count mutations in pruned read
        countMutationsInRead(fields.nucStart, 
                             fields.seqCodesPruned, 
                             refSeq.seqs[target_name],
                             fields.quals, 
                             min_phred,
                             prunedCounts);
        //cout << "mutations counted in pruned seqCodes\n" << flush;
        // count mutations in unpruned read
        countMutationsInRead(fields.nucStart, 
                             fields.seqCodesUnpruned, 
                             refSeq.seqs[target_name], 
                             fields.quals,
                             min_phred,
                             unprunedCounts);
        //cout << "mutations counted in unpruned seqCodes\n" << flush;
    }

    //cout << prunedDiscardedPhredCount << " mutations did not meet phred score filter\n" << flush;

    writeCSV(file_out_path, sample_name, target_name, refSeq, unprunedCounts, prunedCounts);

    return 0;
}
