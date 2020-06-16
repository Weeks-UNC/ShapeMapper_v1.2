// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// Parses SAM alignment file into simpler mutation strings.
// Identifies and removes ambiguously-aligned deletions.
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
#include <cstring>
#include <algorithm>
#include <iterator>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include "string_funcs.h"
#include "ref_seq.h"
#include "parsed_reads.h"

using namespace std;

bool combine_strands;
bool deletion_masking;
bool remove_ambig_del;
bool randomly_primed;
bool trim_both_ends;
int primer_length;
int min_map_qual;

enum SamFields
{
    QNAME,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
};

class Strand
{
// for storing and parsing single strands from SAM sequence reads
// (i.e. the R1 or R2 mate-pair from a single read)
public:
    Strand();
    Strand(const string& line, map<string,string>& refSeqs);
    void init(const string& line, map<string,string>& refSeqs);
    int parseStrand();
    string rawLine;
    vector <string> splitLine; // temp, make private?
    string clusterName;
    string targetName;
    int startIndex; // (1-based)
    int endIndex; // (1-based)
    int mappingQuality;
    string cigarString;
    string rawRead;
    string rawQuals;

    string* refSeq;
    int refIndex; // temp, make private (0-based)
    int readIndex; // temp, make private
    int cigarIndex; // temp, make private
    char lastQual; // temp, make private
    vector<int> cigarNums;
    vector<char> cigarOps;
    vector<char> refAlignedSeqCodes;
    vector<char> refAlignedQuals;
    int splitCIGAR();
    int parseCIGAR(); // calls these helper functions vvv
    // ---------------------------
    int alignmentMatch();
    int insertion();
    int deletion();
    int softClipped();
    // ---------------------------

};


Strand::Strand()
{
    map<string,string> blank;
    init("",blank);
}

Strand::Strand(const string& line, map<string,string>& refSeqs)
{
    init(line, refSeqs);
}


void Strand::init(const string& line, map<string,string>& refSeqs)
{
    rawLine = line;
    trimWhitespace(rawLine);
    splitLineTab(rawLine, splitLine);
    int min_fields = 11;
    if (splitLine.size() >= min_fields){
        clusterName = string(splitLine[QNAME]);
        targetName = string(splitLine[RNAME]);
        startIndex = atoi(splitLine[POS].c_str());
        mappingQuality = atoi(splitLine[MAPQ].c_str());
        cigarString = string(splitLine[CIGAR]);
        rawRead = string(splitLine[SEQ]);
        rawQuals = string(splitLine[QUAL]);
        refSeq = &refSeqs[targetName];
    }else{
        clusterName = "";
        targetName = "";
        refSeq = NULL;
    }
}

int Strand::splitCIGAR()
{
    string opChars = "MIDNSHPX=";
    const char* cCigarString = cigarString.c_str();
    do
    {
        const char *start = cCigarString;
        while(*cCigarString && (opChars.find(*cCigarString)==string::npos))
            cCigarString++;
        if(*cCigarString){
            cigarNums.push_back( atoi(string(start,cCigarString).c_str()) );
            cigarOps.push_back(char(*cCigarString));
        }
    } while (0 != *cCigarString++);
    //cout << "cigarNums length: " << cigarNums.size() << endl << flush;
    //copy(cigarNums.begin(), cigarNums.end(), ostream_iterator<int>(cout, "\t"));
    //cout << "\n" << flush;
    //cout << "cigarOps length: " << cigarOps.size() << endl << flush;
    //copy(cigarOps.begin(), cigarOps.end(), ostream_iterator<char>(cout, "\t"));
    //cout << "\n" << flush;
    return 0;
}


// could maybe do things faster incrementing pointers to char arrays
// or vectors of char
int Strand::alignmentMatch()
{
    //cout << "parsing match\n" << flush;
    int matchLength = cigarNums[cigarIndex];
    int stopReadIndex = readIndex+matchLength;
    char nuc;
    char qual;
    //for i in range(refIndex,refIndex+len(matchLength))
    for (; readIndex<stopReadIndex; readIndex++)
    {
        nuc = rawRead[readIndex];
        qual = rawQuals[readIndex];
        lastQual = qual;
        refAlignedQuals.push_back(qual); // should also be at current refIndex - mapping start offset
        if (nuc!=(*refSeq)[refIndex]){ refAlignedSeqCodes.push_back(nuc); }
        else{ refAlignedSeqCodes.push_back('|'); }
        refIndex++;
    }
    return 0;
}

int Strand::insertion()
{
    //cout << "parsing insert\n" << flush;
    // ignoring for now, maybe store and output to intermediate reads later
    int insertLength = cigarNums[cigarIndex];
    readIndex += insertLength;
    return 0;
}


int Strand::deletion()
{
    //cout << "parsing delete\n" << flush;
    int deletionLength = cigarNums[cigarIndex];
    int i;
    for (i=0; i<deletionLength; i++)
    {
        refAlignedSeqCodes.push_back('-');
        refAlignedQuals.push_back(lastQual); // don't have a phred score within deletion
                                         // (no basecall associated with position), so
                                         // just use last nearby score to approximate
    }
    refIndex += deletionLength;
    return 0;
}

int Strand::softClipped()
{
    //cout << "parsing clipped, cigarIndex = " << cigarIndex << endl << flush;
    int clippedLength = cigarNums[cigarIndex];
    readIndex += clippedLength;
    int i;
    for (i=0; i<clippedLength; i++)
    {
        if (i>=0 && i<(*refSeq).length()){
            refAlignedSeqCodes.push_back('s');
            refAlignedQuals.push_back('#');
            if (cigarIndex == 0) { // leftmost end of read
                startIndex--;
                //cout << "decrementing startIndex\n" << flush;
            }
        }
    }
    return 0;
}

int Strand::parseCIGAR()
{
    if (refSeq){
        //cout << "in parseCIGAR, string(*refSeq) is: \n" << string(*refSeq) << endl << flush;
        //cout << "in parseCIGAR, rawRead is: \n" << rawRead << endl << flush;
        refIndex = startIndex-1; // start parsing at the leftmost matching position (0-based)
        readIndex = 0;
        cigarIndex = 0;
        // seqCodes will store read information relative to each reference nuc
        // ~ = not covered or obscured by long deletion, | = match, ATGC = mismatch,
        // - = 3-prime end of deletion
        while(cigarIndex < cigarOps.size())
        {
            switch(cigarOps[cigarIndex]){
            // codes marked @@@ are not used by bowtie v2.1, probably should edit to support anyway
            case 'M': Strand::alignmentMatch(); break;                
            case '=': break; // @@@
            case 'X': break;// @@@
            case 'I': Strand::insertion(); break;
            case 'D': Strand::deletion(); break;
            case 'N': break; // @@@
            case 'H': break; // @@@
            case 'P': break; // @@@
            case 'S': softClipped(); break;
            default: break; // should raise an exception if we get here
            }
            cigarIndex++;
            //cigarOps.erase(cigarOps.begin());
            //cigarNums.erase(cigarNums.begin());
    }
    }
    endIndex = startIndex + refAlignedSeqCodes.size() - 1;
    //cout << "in parseCIGAR, seqCodes are now: \n"<< string(refAlignedSeqCodes.begin(), refAlignedSeqCodes.end()) << endl << flush;
    return 0;
}

int Strand::parseStrand()
{
    // split CIGAR string
    splitCIGAR();
    // parse CIGAR string and map each basecall to reference sequence
    parseCIGAR();
    return 0;
}


class Read
{
public:
    Read();
    Read(const Strand& R1, const Strand& R2);
    string* refSeq;
    string clusterName;
    int startIndex; // leftmost nuc in reference seq that read maps to (1-based)
    int endIndex;
    vector<char> refAlignedSeqCodes;
    vector<char> refAlignedQuals;
    vector<int> overlapInfo;
    int pruneMutations();
    int removePrimerRegion(int primerLength);
    int testChangeSeqCodes();
};

Read::Read(){}

Read::Read(const Strand& R1, const Strand& R2)
{
    // make a new strand by combining two strands
    // favor higher-quality basecalls when there is a discrepancy and quality difference,
    // otherwise favor the first strand
    //cout << "combining strands\n" << flush;
    refSeq = R1.refSeq;
    clusterName = R1.clusterName;
    //cout << "cluster " << clusterName << "\n" << flush;
    if (!R2.refSeq){ // no R2, just copy important parts of R1
        //cout << "no R2, copying R1\n" << flush;
        startIndex = R1.startIndex;
        endIndex = R1.endIndex;
        refAlignedSeqCodes = R1.refAlignedSeqCodes;
        refAlignedQuals = R1.refAlignedQuals;
        overlapInfo.resize(refAlignedSeqCodes.size());
        int i;
        for (i=R1.startIndex; i<R1.endIndex; i++)
        {
            char c = R1.refAlignedSeqCodes[i-startIndex];
            if(c!='s'){ overlapInfo[i-startIndex] = 1; }
            else{ overlapInfo[i-startIndex] = 0; }
        }
        //fill(overlapInfo.begin(),overlapInfo.end(),1);
    }else{
        //cout << "combining R1 and R2\n" << flush;
        startIndex = min(R1.startIndex, R2.startIndex);
        //cout << "startIndex: " << startIndex << "\n" << flush;
        endIndex = max(R1.endIndex,R2.endIndex);
        //cout << "endIndex: " << endIndex << "\n" << flush;
        if (R1.endIndex<R2.startIndex){
	    //cout << "non-overlapping mate pairs\n" << flush;
        }
        //cout << "len refSeq: "<< (*refSeq).length() << "\n" << flush;
        refAlignedSeqCodes.resize(endIndex-startIndex+1);
        refAlignedQuals.resize(endIndex-startIndex+1);
        overlapInfo.resize(endIndex-startIndex+1);
        int i;
        for (i=startIndex; i<=endIndex; i++) // 1-based
        {
            bool inR1;
            bool inR2;
            inR1 = (i>=R1.startIndex && i<=R1.endIndex);
            inR2 = (i>=R2.startIndex && i<=R2.endIndex);
            char c1;
            char c2;
            char c;
            int q1;
            char q1c;
            int q2;
            char q2c;
            char q;
            int strandCount = 0;
            if(inR1)
            {
                c1=R1.refAlignedSeqCodes[i-R1.startIndex];
                q1c = R1.refAlignedQuals[i-R1.startIndex];
                q1=(int)q1c-33; 
            }
            if(inR2)
            {
                c2=R2.refAlignedSeqCodes[i-R2.startIndex]; 
                q2c = R2.refAlignedQuals[i-R2.startIndex];
                q2=(int)q2c-33;
            } 

            //cout << "got c1, c2, q1, q2\n" <<flush;
            if(inR1 && inR2) // overlap between strands
            {
                //cout << "in overlap region\n" << flush;
                if (c1==c2){ c = c1; q=q1c; strandCount = 2; if(c=='s'){strandCount = 0;}}
                else
                { // strands do not agree
            if (c1=='s'){ c = c2; q=q2c; strandCount = 1;}
                    else if (c2=='s'){ c = c1; q=q1c; strandCount = 1;}
                    else if (q2>q1){ c = c2; q=q2c; strandCount = 1;} // how to count this?
                    else{ c = c1; q=q1c; strandCount = 1;} // fall back to R1 if quals are equal but nucs are different
                }
            }else if(inR1) // only covered by strand 1
            {
                //cout << "in R1 only\n" << flush;
                c = c1;
                q = q1c;
                if (c1!='s'){ strandCount = 1; }
                
            }else if(inR2) // only covered by strand 2
            {
                //cout << "in R2 only\n" << flush;
                c = c2;
                q = q2c;
                if (c2!='s'){ strandCount = 1; }

            }else // between strands, covered by neither
            {
                //cout << "between strands\n" << flush;
            c = '~'; q = '!'; strandCount = 0;
            }
            refAlignedSeqCodes[i-startIndex] = c;
            refAlignedQuals[i-startIndex] = q;
            overlapInfo[i-startIndex] += strandCount;
        }
    }
}

int Read::removePrimerRegion(int primerLength)
{
    // blank out right-most primer region of read (and 1 further nuc) with '~' chars
    int rightZeroIndex = refAlignedSeqCodes.size()-1;
    // to match previous behaviour, if a read overlaps the right end of the reference
    //    sequence, consider the read to start at that nucleotide
    if (endIndex > (*refSeq).length()){ rightZeroIndex -= endIndex-((*refSeq).length()); }
    int leftZeroIndex = max(0,(int)rightZeroIndex - primer_length - 1);
    int i;
    for (i=rightZeroIndex; i > leftZeroIndex; i-- )
    {
        refAlignedSeqCodes[i] = '~';
    }
    // optionally blank out the left end of the read as well
    if (trim_both_ends==true)
    //if (false)
    {
        leftZeroIndex = 0;
        rightZeroIndex = min((int)refAlignedSeqCodes.size()-1,leftZeroIndex + primer_length + 1);
        for (i=leftZeroIndex; i<rightZeroIndex; i++)
	{
	    refAlignedSeqCodes[i] = '~';
        }
    }
    return 0;
}

int Read::pruneMutations()
{
    // replace consecutive mutated nucleotides with a single code at the rightmost (downstream) nuc;
    // identify and remove ambiguously-aligned deletions

    // TODO: break this function into components (ambiguous deletion removal, consecutive mutation replacement, long deletion "masking")
    // Leaving it in its current form for now to maintain consistency with previous version of pipeline.
    vector<char> filteredCodes(refAlignedSeqCodes);

    int totalMismatches = 0;
    int totalDeletions = 0;
    int unambiguousDeletions = 0;
    
    // scan left to right, first nuc to the nuc before last
    int i;
    int k;
    int localIndex;
    char code;
    char downstreamCode;
    for (i=startIndex; i<endIndex; i++)
    {
        localIndex = i-startIndex;
        code = refAlignedSeqCodes[localIndex];
        downstreamCode = refAlignedSeqCodes[localIndex+1];
        if (charFoundIn(code,"ATGC") && !charFoundIn(downstreamCode,"ATGC-")) // if ref nuc code in "ATGC" and downstream nuc code not in "ATGC-"
        { // mismatch with no adjacent downstream mutation
            totalMismatches++;
        }
        else if(code=='-')
        { // deleted nuc
            if(!charFoundIn(downstreamCode,"sATGCN-"))
            { // deleted nuc with no adjacent downstream mutation (right-most deleted nuc)
                totalDeletions++;
                // Ensure deletion is unambiguously aligned
                //  - This isn't guaranteed to work in highly repetitive regions, 
                //    but this takes care of the vast majority of ambiguous deletions. 
                //  - Algorithm amounts to a simple local re-alignment test.   
                
                // find deletion boundaries
                k = 1;
                while(refAlignedSeqCodes[localIndex-k] == '-') k++;
                //cout << "k after finding fivePrimeEnd: " << k << "\n";
                int fivePrimeStart = i-k+1;
                int threePrimeEnd = i;
                //cout << "deletion length " << threePrimeEnd-fivePrimeStart+1 << '\n';
                if (fivePrimeStart-1<0) ;//cout << "fivePrimeStart negative.\n" << flush;
                if (fivePrimeStart-1>=(*refSeq).length()) ;//cout << "fivePrimeStart out of bounds right: "<<fivePrimeStart-1<<"\n" << flush;
                if (threePrimeEnd-1>=(*refSeq).length()) ;//cout << "threePrimeEnd out of bounds right: "<<threePrimeEnd-1<<"\n\n" << flush;
                // slide deletion upstream and downstream as many nucs as the deletion is long
                bool ambiguous = false;
                int offset;
                int offsetRight;
                int offsetLeft;
                string aroundDeletionSeq = (*refSeq).substr(0,fivePrimeStart-1) + (*refSeq).substr(threePrimeEnd, string::npos);
                //cout << "aroundDeletionSeq:\n"<<aroundDeletionSeq << "\n";
                string aroundOffsetDeletionSeq;
                int maxOffset = 1+threePrimeEnd-fivePrimeStart;
                if (remove_ambig_del == true){
                    //bool outOfBounds = false;
                    for (offset= -maxOffset; offset<=maxOffset; offset++)
                    {
                        //cout << "offset " << offset << "\n" << flush;
                        if (offset != 0)
                        {
                            offsetRight = threePrimeEnd+offset+1;
                            offsetLeft = fivePrimeStart+offset;
                            //cout << "offsetRight " << offsetRight <<", offsetLeft " << offsetLeft << "\n" << flush;
                            // if we're over either edge of the reference sequence, ignore this offset
                            if (offsetLeft-1<0  || offsetLeft-1>(*refSeq).length() || 
                                offsetRight-1<0 || offsetRight-1>=(*refSeq).length()){ 
                                    //cout << "offset deletion out of bounds; ignoring\n";
                            }
                            else
                            {
                                aroundOffsetDeletionSeq = (*refSeq).substr(0,offsetLeft-1) + (*refSeq).substr(offsetRight-1, string::npos);
                                //cout << aroundOffsetDeletionSeq << "\n" << flush;
                                // if these sequences match, we've found a simple alternative placement for the deletion
                                if (aroundOffsetDeletionSeq == aroundDeletionSeq)
                                { 
                                    ambiguous = true;
                                    //cout <<"matching sequences around deletion:\n";
                                    //cout << aroundDeletionSeq << "\n";
                                    //cout << aroundOffsetDeletionSeq << "\n" << flush;
                                }
                            }
                        }
                        if(offset==0)
                        {
                            //cout << "deleted sequence: " << (*refSeq).substr(offsetLeft, offsetRight-offsetLeft)<<"\n";
                        }
                    }
                }            

                if (!ambiguous)
                {
                    filteredCodes[localIndex] = '-';
                    unambiguousDeletions++;
                }
                else
                {
                    //cout << "found ambiguous deletion.\n" << flush;
                    // remove this deletion, it is ambiguously-aligned
                    filteredCodes[localIndex] = '|';
                    // remove any deletion-masked nucleotides if ambiguous, since we can't place them accurately
                    offset = -1;
                    while(refAlignedSeqCodes[localIndex+offset] == '-')
                    {
                        filteredCodes[localIndex+offset] = '|';
                        offset--;
                    }
                }
            }
            else
            { // This is a deleted nuc not at the 3 prime end of a long deletion
                if (deletion_masking) {filteredCodes[localIndex] = '~';} // These nucs will not contribute to coverage, since we can't detect a mutation in a sequence we don't have.
                else {filteredCodes[localIndex] = '|';} // Optionally, just treat these nucs as matching the reference, for easier downstream analysis.
            }
        }
    else if (!charFoundIn(code, "~|s"))
    { // This mutation has been pruned, replace with a character indicating a match
        filteredCodes[localIndex] = '|';
    }
    }

    refAlignedSeqCodes.swap(filteredCodes);

    return 0;
}

class SamFile
{
// basic SAM alignment file parsing
// depends on the Strand, Read, and RefSeq classes
public:
    SamFile(const string& filePath);
    Strand R1;
    Strand R2;
    Read combinedRead;
    Read prunedRead;
    Strand Rtemp; // TODO: make this private
    streampos lastStreamIndex; // TODO: make this private
    int nextRead(vector<string>& sequenceTargets, map<string,string>& refSeqs);
    int parseRead();
    string outputRead();
    int printPruningDebug();
    int printDebugInfo();
private:
    ifstream sam;        
    string line;
    int lineCount;
    int nextLine();
    int combineStrands();
    int pruneMutations();
    //int maskLongDeletionCoverage(); // done within pruneMutations() for the moment
    int removePrimerRegion();
};

SamFile::SamFile(const string& filePath)
{
    R1 = Strand();
    R2 = Strand();
    Rtemp = Strand();
    line = "";
    lineCount = 0;
    lastStreamIndex = 0;
    sam.open(filePath.c_str());
    if(sam.is_open()) { //cout << "sam file opened.\n" << flush; 
    }
    else{ 
        stringstream message;
        message << "Unable to open file " << filePath << endl << flush;
        throw runtime_error(message.str().c_str()); 
    }
    
}

int SamFile::nextLine()
{
    // get next non-comment line in sam file and store
    lastStreamIndex = sam.tellg();
    //cout << "lastStreamIndex within nextLine(): " << lastStreamIndex << endl<<flush;
    line = '@';
    while(line[0] == '@' && sam.good())
    {
        //errno = 0;
        getline(sam, line);
        //lineCount++;
        //if (lineCount%100000==0) cout<<"line "<<lineCount<<"\n"<<flush;
    }
    if(!sam.good()){
        /*cerr << "stream not good() after attempting to read line "<<lineCount<<"\n"; 
        cerr << " goodbit: "<< sam.good();
        cerr << " eofbit: "<<sam.eof();
        cerr << " failbit: "<<sam.fail();
        cerr << " badbit: "<<sam.bad()<<"\n"<<flush;
        cerr << "errno: " << strerror(errno) << "\n";
        cerr << "line.length(): "<<line.length() << "\n" <<flush;
        cerr << "lastStreamIndex: "<<lastStreamIndex <<"\n"<<flush;
        sam.clear();
        errno = 0;
        getline(sam, line);
        cerr << " goodbit: "<< sam.good();
        cerr << " eofbit: "<<sam.eof();
        cerr << " failbit: "<<sam.fail();
        cerr << " badbit: "<<sam.bad()<<"\n"<<flush;
        cerr << "errno: " << strerror(errno) << "\n";
        cerr << "line.length(): "<<line.length() << "\n" <<flush;*/
        return 1;
    }
    else { return 0; }
}

int SamFile::nextRead(vector<string>& sequenceTargets, map<string,string>& refSeqs)
{
    // get the next read in the sam file mapping to one of the reference sequences
    // and that is not set as 'N'
    //    - R1 and R2 if combineStrands is set and mapped sequences match
    //    - R1 only otherwise
    R2 = Strand();
    while(true){
        if(nextLine() != 0) { 
            // end of file, stop analysis
            //cout << "end of file, stop analysis\n" << flush;
            return 1;
        }else {
        //cout<<line <<"\n"<< flush;
            Rtemp = Strand(line, refSeqs);
            if (find(sequenceTargets.begin(),
                     sequenceTargets.end(),
                     Rtemp.targetName) != sequenceTargets.end()
                && Rtemp.rawRead != "N"
                && Rtemp.cigarString != "*"
                && Rtemp.mappingQuality >= min_map_qual)
            { 
            //cout << "in nextRead(), Rtemp.rawRead is: " << Rtemp.rawRead << endl << flush;
                //Rtemp.targetName is in reference sequence list and read is not empty
                //cout << "target is in reference sequence list\n" << flush;
                R1 = Rtemp;
                if(combine_strands==true){
                    if(nextLine() != 0){
                        // we've hit the end of file, stop analysis the next time we call this function 
                        //cout << "end of file, stop analysis the next time this function is called\n" << flush;
                        return 0;
                    } 
                    else{ 
                //cout << line<<"\n"<<flush;
                        Rtemp = Strand(line, refSeqs);
                        //cout << "R1.targetName: " << R1.targetName << "\n";
                        //cout << "Rtemp.targetName: " << Rtemp.targetName << "\n";
                        if (   Rtemp.clusterName == R1.clusterName
                            && Rtemp.targetName == R1.targetName 
                            && Rtemp.rawRead != "N" 
                            && Rtemp.mappingQuality >= min_map_qual
			                && abs(Rtemp.startIndex-R1.startIndex)<1000){ // mapped location of R2 is within 1000 nucs of R1
                            R2 = Rtemp;
                            //cout << "R2 target and cluster names matche R1 target and cluster names\n" << flush;
                        }else{ 
                            //seek back to start of line so we can analyze this strand separately the next time we call this function (inefficient but readable)
                            sam.seekg(lastStreamIndex);
                            /*if(!sam.good()){
                                cerr << "stream not good() after attempting to seek\n"; 
                                cerr << " goodbit: "<< sam.good();
                                cerr << " eofbit: "<<sam.eof();
                                cerr << " failbit: "<<sam.fail();
                                cerr << " badbit: "<<sam.bad()<<"\n"<<flush;
                                cerr << "errno: " << strerror(errno) << "\n";
                                exit(0);
                            }*/
                            //cout << "R2 rejected, seeking back to " << lastStreamIndex << endl << flush;
                        }
                    }
                }
                //cout << "got at least one good strand. returning\n" << flush;
                return 0;
            }else{ 
            //cout << "haven't found a strand mapping to a reference sequence yet. continuing\n" << flush;
                continue; // keep looping until we find at least one good strand
            }
        }
    }
}

int SamFile::printPruningDebug()
{
    int offset = min(prunedRead.startIndex,1);
    if (prunedRead.refAlignedSeqCodes.size() > 1){
        cout << "prunedRead.refAlignedSeqCodes:\n" << flush;
        cout << blankString(prunedRead.startIndex-offset) << flush;
        cout << string(prunedRead.refAlignedSeqCodes.begin(), prunedRead.refAlignedSeqCodes.end()) << endl << flush;
    }
    return 0;
}

int SamFile::printDebugInfo()
{
    int offset = min(combinedRead.startIndex,1);

    cout << "R1 start, end: " << R1.startIndex << ", " << R1.endIndex << endl << flush;
    if (R1.refAlignedSeqCodes.size() > 1){
        cout << "R2 start, end: " << R2.startIndex << ", " << R2.endIndex << endl << flush;
    }
    cout << "combined start, end: " << combinedRead.startIndex << ", " << combinedRead.endIndex << endl << flush;
    

    cout << "refSeq:\n";
    cout << blankString(1-offset);
    cout << *(combinedRead.refSeq) << endl << flush;

    if (R1.refAlignedSeqCodes.size() > 1)
    {
        cout << blankString(R1.startIndex-offset) << flush;
        cout << string(R1.refAlignedSeqCodes.begin(), R1.refAlignedSeqCodes.end()) << endl << flush;
        cout << blankString(R1.startIndex-offset) << flush; 
        cout << string(R1.refAlignedQuals.begin(), R1.refAlignedQuals.end()) << endl << flush;
    }

    if (R2.refAlignedSeqCodes.size() > 1)
    {
        cout << blankString(R2.startIndex-offset) << flush;        
        cout << string(R2.refAlignedSeqCodes.begin(), R2.refAlignedSeqCodes.end()) << endl << flush;
        cout << blankString(R2.startIndex-offset) << flush; 
        cout << string(R2.refAlignedQuals.begin(), R2.refAlignedQuals.end()) << endl << flush;
    }    

    cout << "combined seqCodes:\n"; 
    cout << blankString(combinedRead.startIndex-offset) << flush;
    cout << string(combinedRead.refAlignedSeqCodes.begin(), combinedRead.refAlignedSeqCodes.end()) << endl << flush;

    cout << "overlap info:\n";    
    cout << blankString(combinedRead.startIndex-offset) << flush;
    copy(combinedRead.overlapInfo.begin(), combinedRead.overlapInfo.end(), ostream_iterator<int>(cout, ""));
    cout << "\n" << flush;    

    return 0;
}


int SamFile::combineStrands()
{
    combinedRead = Read(R1,R2); 
    //printDebugInfo();    
    return 0;
}

int SamFile::pruneMutations()
{
    prunedRead = Read(combinedRead);
    prunedRead.pruneMutations();
    //printPruningDebug();
    //printDebugInfo();
    return 0;
}

//int SamFile::maskLongDeletionCoverage(){return 0;}
int SamFile::removePrimerRegion()
{
    // TODO: have primer clipping option for each target RNA, instead of current "global" option
    // This function attempts to blank out the right-most reverse transcription primer region of this read with '~' chars
    combinedRead.removePrimerRegion(primer_length);
    return 0;
}

int SamFile::parseRead()
{
    R1.parseStrand();
    R2.parseStrand();
    combineStrands();
    if (randomly_primed) removePrimerRegion();
    pruneMutations();
    //maskLongDeletionCoverage(); // done within pruneMutations() for now    
    return 0;
}

string SamFile::outputRead()
{
    // line format: 
    // startIndex (1-based, may be negative) TAB endIndex (1-based, may be larger than sequence target) TAB prunedSeqCodes TAB un-prunedSeqCodes
    string outputLine;
    outputLine = combinedRead.clusterName + "\t" + 
                 intToString(prunedRead.startIndex) + "\t" +
                 intToString(prunedRead.endIndex) + "\t" +
                 string(prunedRead.refAlignedSeqCodes.begin(), prunedRead.refAlignedSeqCodes.end()) + "\t" +
                 string(combinedRead.refAlignedSeqCodes.begin(), combinedRead.refAlignedSeqCodes.end()) +"\t" +
                 string(combinedRead.refAlignedQuals.begin(), combinedRead.refAlignedQuals.end()) +"\n";  
    //cout << "SamFile::outputRead() called\n" << flush;     
    return outputLine; 
}

static void printUsage()
{
    cerr<<"Usage example:\n    parseAlignment -combine_strands -deletion_masking "
    <<"-randomly_primed -trim_both_ends -primer_length 10 -ref_seqs ref.fa "
    <<"-file_in test.sam -out_folder simplified_mutation_strings\n"<<flush;
}

int main(int argc, char* argv[])
{
    // commandline argument handling --------------------------
    combine_strands = false; // corresponds to the alignPaired config option
    deletion_masking = false;
    randomly_primed = false; // if true, exclude primer region from counting
    trim_both_ends = false; // if we don't know which strand the RNA came from (e.g. transriptome),
                               // we need to remove the RT primer from both ends of each read to be safe
    primer_length = 10;
    min_map_qual = 10;
    string ref_seqs_path;
    string file_in_path;
    string out_folder;
    
    remove_ambig_del = false;

    char* endptr;
    int base = 10;    
    int optionIndex = 0;
    static struct option longOptions[] = {
    {"combine_strands",      no_argument, NULL, 'c' },
    {"deletion_masking",     no_argument, NULL, 'm' },
    {"randomly_primed", no_argument, NULL, 'r' },
    {"trim_both_ends", no_argument, NULL, 'b' },
    {"primer_length",    required_argument, NULL, 'p' },
    {"min_map_qual", required_argument, NULL, 'q' },
    {"ref_seqs", required_argument, NULL, 's' },
    {"file_in", required_argument, NULL, 'i' },
    {"out_folder", required_argument, NULL, 'o' },
    {"remove_ambig_del", no_argument, NULL, 'a' },
    {NULL,    no_argument,            NULL,  0 }
    };
    static const char *optString = ":c:m:r:b:p:q:s:i:o:a";

    static int flags[3] = { 0 };

    while (getopt_long_only(argc, argv, optString, longOptions, &optionIndex) != -1)
    {
        switch(optionIndex)
        {
        case 0:
            errno = 0;
            combine_strands = true;
            break;
        case 1:
            errno = 0;
            deletion_masking = true;
            break;
        case 2:
            errno = 0;
            randomly_primed = true;
            break;
        case 3:
	    errno = 0;
	    trim_both_ends = true;
            break;
        case 4:
            errno = 0;
            primer_length = strtol(optarg, &endptr, base);
            if(errno > 0){ cerr<<"Unable to parse primer_length argument.\n"<<flush; exit(1); }
            break;
        case 5:
            errno = 0;
            min_map_qual = strtol(optarg, &endptr, base);
            if(errno > 0){ cerr<<"Unable to parse min_map_qual argument.\n"<<flush; exit(1); }
            break;
        case 6:
            flags[0] = 1;
            ref_seqs_path = string(optarg);
            break;
        case 7:
            flags[1] = 1;
            file_in_path = string(optarg);
            break;
        case 8:
            flags[2] = 1;
            out_folder = string(optarg);
            if(out_folder[out_folder.length()-1]!='/'){ out_folder = out_folder+"/"; }
            break;
        case 9:
	    errno = 0;
            remove_ambig_del = true;
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
    for (i=0; i<3; i++)
    {
        if (!flags[i]){ foundAllFlags = false; }
    }
    if (!foundAllFlags){
        cerr << "ref_seqs, file_in, and out_folder are required arguments.\n"<<flush;
        printUsage();
        exit(1);
    }
    // --------------------------------------------------------
    //cout << "started.\n" << flush;
    RefSeq refSeq(ref_seqs_path);
    //refSeq.printSequences();

    //cout << "ref seqs loaded.\n" << flush;
    SamFile samFile(file_in_path);

    string sample_name = file_in_path;
    // remove folders and extension
    size_t charIndex = sample_name.rfind('/');
    if (charIndex != string::npos){ sample_name.erase(0, charIndex+1); }
    charIndex = sample_name.rfind('.');
    if (charIndex != string::npos){ sample_name.erase(charIndex); }

    ParsedReadOutputFiles parsedReadFiles(out_folder,
                                          sample_name,
                                          refSeq.seqNames);
    //cout << "starting read input\n" << flush;
    //int readCount = 0;
    while( samFile.nextRead(refSeq.seqNames, refSeq.seqs) == 0){
        samFile.parseRead();
        //readCount++;
        //if (readCount%100==0){ cerr << "read " << readCount << "\n" << flush; }
        //cout << "writing read string\n";
        parsedReadFiles.writeLine(samFile.R1.targetName, samFile.outputRead());
        //cout << samFile.outputRead();
        //cout << "\n" << flush;
        //cout << samFile.R1.clusterName << "\t" << samFile.R1.targetName << endl << flush;
    }

    return 0;
}
