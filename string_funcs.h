// Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
// String processing functions.
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

#include <string>
#include <map>
#include <vector>

std::string intToString(int i)
{
    char chars[14];
    return std::string(chars, sprintf(chars, "%d", i));
}

bool isOnlyWhitespace(const std::string& s)
{
    return (s.find_first_not_of(" \t\n\r")==std::string::npos);
}

bool charFoundIn(const char& c, const std::string& list)
{
    size_t foundIndex;
    foundIndex = list.find_first_of(c);
    if (foundIndex != std::string::npos) return true;
    else return false;
}

int filterChars(std::string& s, const std::string& keepChars)
{
    // remove all chars not found in keepChars from string s (in place)
    size_t foundIndex;
    foundIndex = s.find_first_of(keepChars);
    if (foundIndex == std::string::npos){ s.clear(); return 1;}
    else { while (true)
    {
        foundIndex = s.find_first_not_of(keepChars);
        if (foundIndex != std::string::npos){ s.erase(foundIndex,1); }
        else break;    
    }}
    return 0;
}

std::string blankString(int length)
{
    std::string blanks = "";
    int i;
    for(i=0;i<length;i++)
    {
        blanks += " ";
    }
    return blanks;
}

int trimWhitespace(std::string& s)
{
    // remove whitespace chars from the start and end of string s (in place)
    size_t foundIndex;
    std::string filterChars = " \t\n\r";
    foundIndex = s.find_first_not_of(filterChars);
    if (foundIndex == std::string::npos){ s.clear(); return 1;}
    else { 
	// trim from the left
        foundIndex = s.find_first_not_of(filterChars);
        if (foundIndex != std::string::npos){ s.erase(0,foundIndex); }
        // trim from the right
        foundIndex = s.find_last_not_of(filterChars);
        if (foundIndex != std::string::npos){ s.erase(foundIndex+1); }    
    }
    return 0;
}

int splitLineTab(const std::string& line, std::vector<std::string>& splitLine)
{
    // Couldn't find strtok in my headers, maybe my build env is broken
    // Using this instead.
    const char* cLine = line.c_str();
    do
    {
        const char *start = cLine;
        while(*cLine && *cLine != '\t')
            cLine++;
        splitLine.push_back(std::string(start, cLine));
    } while (0 != *cLine++);
    return 0;
}
