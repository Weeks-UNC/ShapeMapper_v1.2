// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// Functions for loading fasta sequence files.
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
#include <string>
#include <map>
#include <vector>
#include <stdexcept>


// function to correctly handle different line endings
// (taken from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf)
std::istream& safeGetline(std::istream& is, std::string& t)
{
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();

  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\n':
      return is;
    case '\r':
      if(sb->sgetc() == '\n')
	sb->sbumpc();
      return is;
    case EOF:
      // Also handle the case when the last line has no line ending
      if(t.empty())
	is.setstate(std::ios::eofbit);
      return is;
    default:
      t += (char)c;
    }
  }
}

class RefSeq
{
// storage class for reference sequences and associated names
public:       
    RefSeq(const std::string& filePath);
    std::vector <std::string> seqNames; // will hold copies of keys to seq map
    std::map <std::string, std::string> seqs;
    std::string keepNucs;
    int printSequences();
};

RefSeq::RefSeq(const std::string& filePath)
{
    keepNucs = "ATGC";
    std::string line;
    std::ifstream fa(filePath.c_str());
    if (fa.is_open())
    {
        int seqIndex = -1;
        std::string currentSeqName;
        while(fa.good())
        {
	    getline(fa, line);
	    //std::cout << "ref seq file tellg(): " << fa.tellg() << std::endl << std::flush;
            // test for fasta sequence name
            if (line[0] == '>'){ 
                seqIndex++;
                trimWhitespace(line);
                line.erase(0,1); // remove leading '>' char from seq name
                currentSeqName = line;
                seqNames.push_back(currentSeqName);
                seqs[currentSeqName] = "";
            }else{
                filterChars(line, keepNucs);  
                seqs[currentSeqName].append(line);
            }
	}        
        fa.close();
    }else{ 
        std::stringstream message;
        message << "Unable to open file " << filePath << std::endl << std::flush;
        throw std::runtime_error(message.str().c_str()); 
	//cout << "Unable to open file " << filePath << endl;
    }

}

int RefSeq::printSequences()
{
    int i;
    for (i=0; i<seqNames.size(); i++)
    {
        std::string seqName = seqNames[i];
        std::string seq = seqs[seqName];
        std::cout << "name: " << seqName << "\n" << "seq: " << seq << "\n\n" << std::flush;
    } 
    return 0;
}
