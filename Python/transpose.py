#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2014 Junli Zhang <zhjl86@gmail.com>
# Function: transpose a file (row to col) 
# Default is white space delimited
# Need one argument for filename.
# This program is free software but WITHOUT ANY WARRANTY
# 


import sys

if len(sys.argv) != 2:
	print "Error: You should specify the input file!!!"
	exit()

infile = sys.argv[1]
# read files and split each line to a list
with open(infile) as f:
	lis=[x.split('\t') for x in f] # list of list


out = open("transposed_" + infile, "w")
for x in zip(*lis): # unzip a list
	out.write('\t'.join(map(str,x)) + '\n')
