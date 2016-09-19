#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
#  
# Version 3: modified the script so it can handle more than 3 environments
#            now it put the pvalues of all markers in a dict, then filter before writing.

from glob import glob

# change here for different trait and different NAMs
traits = ["TotalSN", "Yield", "KW400", "SW6spikes", "KW6spikes", "Seeds6spikes", "Spikes.sqft"]
nams = ["RAC875", "RSI5", "UC1036", "UC1419"]

# the basic function
def sigSNP (trait, nam):
	raw = glob("pvalues*" + nam + "*")
	raw.sort()

	out = open(trait + "_" + nam + ".txt", "w")
	out.write("SNP") # the first line of output
	pvalues = {}

	for ff in raw:
		out.write("\t" + ff)
		with open(ff) as infile:
			#next(infile) # skip header
			ID = infile.readline().rstrip() # the first line is the trait names
			n1 = ID.split('\t').index(trait) # the column of the trait
			for line in infile:
				col = line.rstrip().split('\t')
				#key = ",".join(col[0:3])
				key = col[0]
				#if float(col[n1]) > 0.01: continue # skip the line if p> 0.01
				pvalues.setdefault(key, []) # avoid key error for the first time
				pvalues[key].append(col[n1]) # values are a list
				
				

	out.write("\n")
	for k, v in pvalues.iteritems():
		if sum(float(i) < 0.01 for i in v) > 1: # whether at least 2 values < 0.01
			out.write(k + "\t" + "\t".join(v) + "\n")
	out.close()

## loop through traits and NAMs

for x in traits:
	for y in nams:
		sigSNP(x, y)
