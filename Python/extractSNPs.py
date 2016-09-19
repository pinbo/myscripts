#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
#  
from glob import glob

# change here for different trait and different NAMs
trait = "TotalSN"
nam = "DD"

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
			if float(col[n1]) > 0.01: continue # skip the line if p> 0.01
			pvalues.setdefault(key, "") # avoid key error for the first time
			pvalues[key] += '\t' + col[n1]

out.write("\n")
for k, v in pvalues.iteritems():
	if v.count('\t') > 1:
		out.write(k + v + "\n")
out.close()
