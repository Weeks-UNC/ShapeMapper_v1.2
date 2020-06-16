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

import csv,sys,os

def pivot(infile,outdir):

    with open(infile) as f:
        reader = csv.reader(f)
        cols = []
        for row in reader:
            cols.append(row)
    outfile = os.path.join(outdir,os.path.splitext(os.path.split(infile)[1])[0]+".csv")
    with open(outfile, 'wb') as f:
        writer = csv.writer(f)
        for i in range(len(max(cols, key=len))):
            writer.writerow( [(c[i] if i<len(c) else '') for c in cols])
            
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: <inputFile> <outputDir>"
        exit(0)
    pivot(sys.argv[1], sys.argv[2])
