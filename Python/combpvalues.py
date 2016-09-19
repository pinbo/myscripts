#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  combpvalues.py
#  
#  Copyright 2013 Junli Zhang <zhjl86@gmail.com>
#  Extract p value columns and combine them from GAPIT OUTPUTs

from glob import glob
raw = glob("*.GWAS.Results.csv")
raw.sort()

out = open("combpvalues.csv", "w")
trait = [gwas[7:-17] for gwas in raw]
out.write("SNP,Chrom,Position," + ','.join(trait) + "\n")
markers = {}

for gwas in raw:
    with open(gwas) as infile:
        next(infile) # skip header
        for line in infile:
            col = line.split(',')
            key = ",".join(col[0:3])
            markers.setdefault(key, "") # avoid key error for the first time
            markers[key] += "," + col[3]

for k, v in markers.iteritems():
    out.write(k + "," + v + "\n")
out.close()

