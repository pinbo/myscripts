#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  prepare_rqtl.py
#  
#  Copyright 2015 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# Function: combine genotyping data and maps for rqtl analysis

#from glob import glob
import sys

if len(sys.argv) != 4:
	print 'Usage: prepare_rqtl.py Genoptype_file Map_file Outputfile'
	sys.exit()
#print sys.argv[1]
genofile, mapfile = sys.argv[1:3]
out = open(sys.argv[3], "w")

marker = {} # dictionary for all marker postion and genotype
# read map file
with open(mapfile) as infile:
	#next(infile) # skip header
    for line in infile:
		col = line.rstrip().split()
		#print col
		key = col[1].split(":")[0] # some maps have RAC875_c42700_264:1AS
		marker[key] = ",".join([col[0], col[2]])

# read genoptying file
n = 0 # unknown marker position
dd = {} # just to save markers that have genotypes
with open(genofile) as infile:
	next(infile) # skip header
	for line in infile:
		col = line.replace("\t", ",").split(",", 1)
		key = col[0]
		dd[key] = 0 # just to save keys for late use
		if key in marker:
			marker[key] += "," + col[1]
		else:
			marker[key] = ",".join(["unknown", str(n), col[1]])
			n += 0.1

# I found some markers have map positions but do not have genoptyes
# so only to save those with genotype data			
for k in dd:
    out.write(k + "," + marker[k])
out.close()
