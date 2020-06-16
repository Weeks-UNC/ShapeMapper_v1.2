// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// Simple fastq read quality trimming program.
// Copyright Steven Busan 2014

/*---------------------------------------------------------------------------
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
*///--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h> // character handling functions
#include <getopt.h> // cmdline long arg parsing

// these variables will be overwritten by commandline args passed to this program
int MINPHRED = 10;
int MINLENGTH = 25;
int WINDOWSIZE = 4;
int MAXLINELENGTH = 5000; 

static int printEntireArray(unsigned char *a)
{
    // for debugging purposes, print an entire char array, including null chars
    int i;
    unsigned char c;
    for(i = 0; i < MAXLINELENGTH; i++)
    {
        c = a[i];
        if (c=='\0'){ c='_'; } // replace null chars with underscores
	printf("%c", c);
    }
    printf("\n\n");
    return 0;
}

static size_t getFileSize(const char *filePath)
{
    struct stat sb;
    if (stat(filePath, &sb) != 0)
    {
        fprintf(stderr, "Failed to read file %s: %s", filePath, strerror(errno));
       	exit(1);
    }
    //fprintf(stdout, " sizeof(sb.st_size) = %zu\n", sizeof(sb.st_size));
    //fprintf(stdout, " file size = %zu\n", sb.st_size); // this value is correct
    //size_t fileSize = sb.st_size;
    return sb.st_size;
}

typedef struct
{
    // output file handle and path
    const char* filePath;
    FILE* f;
}
OutputFile;

static void openOutputFile(OutputFile* a)
{
    a->f = fopen(a->filePath,"w");
    if (!a->f)
    {
	fprintf(stderr,"Unable to read file %s: %s", a->filePath, strerror(errno));
	exit(1);
    }
}

static void closeOutputFile(OutputFile* a)
{
    // close file f.
    int status;
    status = fclose(a->f);
    if (status!=0)
    {
        fprintf(stderr, "Error closing file %s: %s", a->filePath, strerror(errno));
	exit(1);
    }
}

typedef struct
{
    // input file handle, path, and associated memory block for buffered read
    const char* filePath;
    FILE* f;
    size_t fileSize;
    size_t blockSize;
    size_t blockIndex;
    char* blockContents;
} 
InputFile;

static void openInputFile(InputFile* a)
{
    // open the file given by filePath, return a handle to it.
    // also allocate a block of memory for buffered read.
    // set starting block index to 0.
    
    //printf("filePath: %s\n",a->filePath);

    // get file size so we know how much memory to allocate
    a->fileSize = (size_t)getFileSize(a->filePath)+(size_t)1; // +1 leaves room for null char apparently
    //printf("fileSize in openInputFile(): %zu\n",a->fileSize);
    //printf("size returned from getFileSize: %zu\n", getFileSize(a->filePath));

    a->blockSize = (size_t)4096;
    if (a->fileSize < a->blockSize){ a->blockSize=a->fileSize; }

    a->blockContents = malloc(a->blockSize);

    if (!a->blockContents)
    {
	fprintf(stderr, "Unable to read file %s: Out of memory.", a->filePath);
        exit(1);
    }
    a->f = fopen(a->filePath,"r");
    if (!a->f)
    {
	fprintf(stderr,"Unable to read file %s: %s", a->filePath, strerror(errno));
	exit(1);
    }
    a->blockIndex = 0;
}

static void readNextBlock(InputFile* a)
{
    // Alternatively, could pipe file to stdin, since stdin is already buffered.
    // This would require one more process to stay active to do the reading from disk
    size_t bytesReturned;
    bytesReturned = fread(a->blockContents, sizeof(unsigned char), a->blockSize, a->f);
    //printf("blockIndex %zu, blockSize %zu, fileSize %zu\n",a->blockIndex,a->blockSize,a->fileSize);
    size_t expectedBytes = a->blockSize; 
    if (a->fileSize < a->blockSize*(a->blockIndex+(size_t)1)){ // this file size is not always correct
        expectedBytes = a->blockSize - (a->blockSize*(a->blockIndex+(size_t)1) - a->fileSize) - (size_t)1;
    }
    //printf("%zu bytes expected\n",expectedBytes);
    //printf("%zu bytes returned from fread\n\n",bytesReturned);
    
    if (bytesReturned != expectedBytes)
    { 
        fprintf(stderr, "Error reading file %s: %s\n", a->filePath, strerror(errno));
        fprintf(stderr, " %zu bytes expected, got %zu\n", expectedBytes, bytesReturned);
        fprintf(stderr, " blockSize %zu, blockIndex %zu, fileSize %zu\n", a->blockSize, a->blockIndex, a->fileSize);
        fprintf(stdout, " sizeof(fileSize) = %zu", sizeof(a->fileSize));
        fprintf(stderr, a->blockContents);
	exit(1);
    }
    //printf("%s",a->blockContents);
    a->blockIndex++;    
}

static void closeInputFile(InputFile* a)
{
    // close file f.
    // free memory block associated with buffered read.
    int status;
    status = fclose(a->f);
    if (status!=0)
    {
        fprintf(stderr, "Error closing file %s: %s", a->filePath, strerror(errno));
	exit(1);
    }
    free(a->blockContents);
}


static int charToPhred(char c)
{
    int asciiValue = (int)c;
    int phredScore = asciiValue - 33;
    return phredScore;
}

static int locateFirstBadNuc(char *qualities)
{
    int badNucIndex = -1;
    unsigned int i=0;
    char c = ' ';
    while (i<MAXLINELENGTH && c != '\n')
    {
        c = qualities[i];  // could try pointer increment and see if faster?
        if (charToPhred(c) < MINPHRED){
            badNucIndex = i;
            break;
        }
        i++;
    }
    return badNucIndex; // >= 0 indicates bad nuc index, -1 indicates no bad nuc found
}

static int locateFirstBadNucWindowed(char *qualities)
{
    unsigned int scores[MAXLINELENGTH];
    unsigned int i=0;
    char c = ' ';
    while (i<MAXLINELENGTH && c != '\n')
    {
        c = qualities[i];
        scores[i] = charToPhred(c);
        i++;
    }

    int sum = 0;
    int mean = 0;
    int badNucIndex = -1;
    i=0;
    c = ' ';
    unsigned int j=0;
    while (i<MAXLINELENGTH-WINDOWSIZE && c != '\n')
    {
        sum=0;
        for (j=0; j<WINDOWSIZE; j++){
	  sum += scores[i+j];
        }
        mean = sum/WINDOWSIZE;
        if (mean < MINPHRED){
            badNucIndex = i;
            break;
        }
        i++;
    }
    return badNucIndex; // >= 0 indicates bad nuc index, -1 indicates no bad nuc found
}

const char *byteToBinary(int x)
{
    static char b[9];
    b[0] = '\0';

    int z;
    for (z = 128; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}

static int filterFileContents(InputFile* a, OutputFile* b)
{
    char lineBlock[4][MAXLINELENGTH];

    // scan and store each char in appropriate temp line buffer until newline character
    size_t localLineIndex = 0; // 0-3, line index relative to current sequencing read
    size_t columnIndex = 0;
    size_t readIndex = 0; // index of current read, 0-based
    size_t charIndex = 0; // index of current char in file contents buffer
    size_t localIndex = 0; // index of current char within block (when >= blocksize, we need a new block)
    char currentChar = ' ';
    int badNucIndex;

    // get the first block from the file
    readNextBlock(a);

    //printf("charIndex %i, *inputSize %zu\n", charIndex, *inputSize);
    while (charIndex < a->fileSize-1 && 
           currentChar != EOF)
    {
        // loop until newline or EOF
        while (currentChar != '\n' && currentChar != EOF && charIndex < a->fileSize-1)
	{	
            currentChar = a->blockContents[localIndex];
	    charIndex++;
            localIndex++;
	    lineBlock[localLineIndex][columnIndex] = currentChar;
            columnIndex++;

            // get a new block from the file if we need it
            if (localIndex >= a->blockSize)
	    {
                //printf("charIndex %i\n",charIndex);
                readNextBlock(a);
                localIndex = 0;
            }
	}
        lineBlock[localLineIndex][columnIndex] = '\0'; // terminate string
        currentChar = ' ';
        localLineIndex++;
        columnIndex = 0;
        if (localLineIndex > 3)
	{
            //printf("read %u, char %u\n",readIndex, charIndex);
            localLineIndex = 0;
            readIndex += 1;
            //printf("%s%s\n\n",lineBlock[1],lineBlock[3]);            
            // find a nuc (if any) to trim past
            //badNucIndex = locateFirstBadNuc(lineBlock[3]);
            badNucIndex = locateFirstBadNucWindowed(lineBlock[3]);
            if (badNucIndex == -1){}
            else if (badNucIndex < MINLENGTH){
	    // trim entire read if remaining read does not pass the length cutoff
		lineBlock[1][0] = 'N';
		lineBlock[1][1] = '\n';
                lineBlock[1][2] = '\0'; // terminate string
		lineBlock[3][0] = '#';
		lineBlock[3][1] = '\n';
                lineBlock[3][2] = '\0'; // terminate string
	    }else{
            // otherwise trim read downstream of the first bad basecall
		lineBlock[1][badNucIndex] = '\n';
                lineBlock[1][badNucIndex+1] = '\0';
                lineBlock[3][badNucIndex] = '\n';
                lineBlock[3][badNucIndex+1] = '\0';
	    }
            //printf("%s%s%s%s",lineBlock[0],lineBlock[1],lineBlock[2],lineBlock[3]);            
            //printEntireArray(lineBlock[1]);
            int i;
            for(i=0; i<4; i++) // loop over 4 lines in read (name, sequence, delimiter, qualities)
            { 
                if (fputs(lineBlock[i],b->f) < 0 ){ // write line to output file
                    fprintf(stderr, "Error writing to file %s: %s", b->filePath, strerror(errno));
	            exit(1);
                }
            }
        }
           
    }
    //printf("%u total chars read\n",charIndex+1);

    return 0;
}

static void printUsage()
{
    fprintf(stderr,"Usage example (all args required):\n    filterPhred -minphred 10 -minlength 25 -maxlinelength 5000 -filein test.fastq -fileout filtered_reads/test.fastq\n");
}

int main(int argc, char* argv[])
{
    printf("Starting filtering\n");

    InputFile fileIn;
    OutputFile fileOut;
    char* endptr;
    int base = 10;    

    int optionIndex = 0;
    static struct option longOptions[] = {
	{"minphred",      required_argument, NULL, 'p' },
        {"windowsize",    required_argument, NULL, 'w' },
	{"minlength",     required_argument, NULL, 'l' },
	{"maxlinelength", required_argument, NULL, 'b' },
	{"filein",    required_argument, NULL, 'i' },
	{"fileout",   required_argument, NULL, 'o'},
	{NULL,    no_argument,            NULL,  0 }
    };
    static const char *optString = ":p:l:b:i:o:h";

    static int flags[6] = { 0 };
    int lastOptionIndex = -1;
    while (getopt_long_only(argc, argv, optString, longOptions, &optionIndex) != -1){
        if(optionIndex==lastOptionIndex){ fprintf(stderr,"Error reading commandline arguments.\n"); exit(1); printUsage();}
        lastOptionIndex = optionIndex; // TODO: generate a more helpful error message when an unrecognized arg is found
                                       // haven't been able to figure out why getopt_long_only is behaving as it is
        switch(optionIndex){
        case 0:
            flags[0] = 1;
            errno = 0;
            MINPHRED = strtol(optarg, &endptr, base);
            if(errno > 0){ fprintf(stderr,"Unable to parse minPhred argument.\n"); exit(1); }
            break;
        case 1:
            flags[1] = 1;
            errno = 0;
            WINDOWSIZE = strtol(optarg, &endptr, base);
            if(errno > 0){ fprintf(stderr,"Unable to parse windowSize argument.\n"); exit(1); }
            break;
        case 2:
            flags[2] = 1;
            errno = 0;
            MINLENGTH = strtol(optarg, &endptr, base);
            if(errno > 0){ fprintf(stderr,"Unable to parse minLength argument.\n"); exit(1); }
            break;
        case 3:
            flags[3] = 1;
            errno = 0;
            MAXLINELENGTH = strtol(optarg, &endptr, base);
            if(errno > 0){ fprintf(stderr,"Unable to parse maxLineLength argument.\n"); exit(1); }
            break;
        case 4:
            flags[4] = 1;
            fileIn.filePath = optarg;
            break;
        case 5:
            flags[5] = 1;
            fileOut.filePath = optarg;
            break;
        case 'h':
        case '?':
	    // getopt_long_only() doesn't seem to ever set this value
	    fprintf(stderr,"Unrecognized argument \"%s\"\n",optopt);
        default:
            printUsage();
            exit(0);
            break;
        }
    }
    int foundAllFlags = 1; // True
    size_t i;
    for (i=0; i<6; i++)
    {
        //printf("arg %zu flag is %i\n",i,flags[i]);
        if (!flags[i]){ foundAllFlags = 0; }
    }
    if (!foundAllFlags){
        fprintf(stderr,"All arguments are required.\n");
        printUsage();
        exit(1);
    }

    openInputFile( &fileIn );
    openOutputFile( &fileOut );

    filterFileContents( &fileIn, &fileOut );

    closeInputFile( &fileIn );
    closeOutputFile( &fileOut );
    printf("    Filtered %s and wrote output to %s\n",fileIn.filePath, fileOut.filePath);
    printf("Filtering complete\n");
    return(0);
}
