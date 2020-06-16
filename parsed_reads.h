// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// Functions for handling mutation string files.
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

class ParsedReadOutputFiles
{
// stores handles to files for outputing parsed sequence mutation codes
public:
    ParsedReadOutputFiles(const std::string& outputFolder, 
                          const std::string& sampleName, 
                          const std::vector<std::string>& sequenceTargets);
    std::map<std::string,std::ofstream*> parsedReadFiles;
    int writeLine(const std::string& targetName, const std::string& line);
};

ParsedReadOutputFiles::ParsedReadOutputFiles(const std::string& outputFolder, 
                                             const std::string& sampleName, 
                                             const std::vector<std::string>& targetNames)
{
    std::string filePath;
    int i;
    for (i=0; i<targetNames.size(); i++){
        filePath = outputFolder+sampleName+"_"+targetNames[i]+".txt";
        parsedReadFiles[targetNames[i]] = new std::ofstream(filePath.c_str());
        if (!(*(parsedReadFiles[targetNames[i]])).is_open())
        {
	    std::stringstream message;
	    message << "Unable to open file " << filePath << std::endl << std::flush;
	    throw std::runtime_error(message.str().c_str()); 
        }
    }
}

int ParsedReadOutputFiles::writeLine(const std::string& targetName, const std::string& line)
{
    *(parsedReadFiles[targetName]) << line << std::flush;
    // Following unformatted output may be faster, but seems to silently fail to output
    // a small number of lines for some reason.
    //(*(parsedReadFiles[targetName])).write( line.c_str(), sizeof(char)*line.size() );
    return 0;
}


    
