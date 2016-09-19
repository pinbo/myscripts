#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  combine_indices_based_plot.py
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

from glob import glob

env = "Davis2015" # change here for different environments
raw = glob("Indices*") # All file names start from "Indices"
raw.sort()

out = open("condensed_indices_" + env + ".txt", "w")
#out.write("Wavelength," + ','.join(raw) + "\n")
csr = {}

for ff in raw:
    with open(ff) as infile:
        #next(infile) # skip header
        for line in infile:
			# remove double quotes first, then strip newline, then split once
            col = line.replace('"', '').rstrip().split('\t', 1)
            # I can also use: key, value = line.split()
            key = col[0]
            csr.setdefault(key, "") # avoid key error
            csr[key] += "\t" + col[1]

# sort characters before numbers
# only for integers: if it is an integar then change it to integars
# otherwise change it to 0, so it is on top after sorting
def f(e):
	if e.isdigit():
		return int(e) or 0 

# for complex numbers
#def is_number(s):
    #try:
        #float(s)
        #return True
    #except ValueError:
        #return False

#for k in sorted(csr, key = lambda x: int(x)):
for k in sorted(csr, key=f):
    out.write(k + csr[k] + "\n")
out.close()
