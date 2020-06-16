# Part of the SHAPE-MaP data analysis pipeline (ShapeMapper).
# Loads stored configuration options.
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

import os, pickle

# assuming parseConfig.py has already been run, 
# load the pickled dictionary of parameters
fileIn = open(os.path.join(os.getcwd(),"temp_config.pickle"), "r")
unpickled = pickle.load(fileIn)
fileIn.close()

#print str(unpickled)

# move the conf dictionary into the global namespace
for key in unpickled.keys():
    globals()[key] = unpickled[key]
del unpickled
