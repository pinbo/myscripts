#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  Copyright 2020 Junli Zhang <zhjl86@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.

# convert Windows cartographer eQTL file to organized csv file
# two arguments: 1. input eqtl output file from Windows Cartographer (file name is like xxxx_out_in-c-eqtl.txt), 2. output.csv

# At the QTL results window, run "automatic locating QTL", then if you can save result as Excel with two tabs, then you do not need this script
# I have a problem saving the QTL results into the Excel file, so I save it as EQTL file, then using this script to organize the result

#from glob import glob
import sys
import re

if len(sys.argv) != 3:
	print 'Usage: ./convert-cartographer-qtl.py organized-result.csv'
	sys.exit()
#print sys.argv[1]
qtlfile = sys.argv[1]
out = open(sys.argv[2], "w")
out.write("Trait,number,Chrom,Markr,Position,LOD,Additive,Dominance,R2,TR2,S\n")
# read qtl file
trait = "xx"
with open(qtlfile) as infile:
	#next(infile) # skip header
    for line in infile:
		if line.startswith("# for trait"):
			trait = line.rstrip().split()[4]
		elif line.startswith("#") or line.startswith("-"):
			continue
		else:
			out.write(re.sub(' +', ',', trait + line))

out.close()