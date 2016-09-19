#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  condense_reference.py
#  Function: combine dark and white references for adjustment
#
#  Copyright 2014 Junli Zhang <zhjl86@gmail.com>
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
import sys

channel = sys.argv[1]

raw = glob(channel + "*") # "Sdark0" for example
raw.sort()

out = open("condensed_" + channel + ".csv", "w")
out.write("Wavelength," + ','.join(raw) + "\n")
csr = {}

for reference in raw:
    with open(reference) as infile:
        #next(infile) # skip header
        for line in infile:
            col = line.split() # spaces
            # I can also use: key, value = line.split()
            key = col[0]
            csr.setdefault(key, "") # avoid key error
            csr[key] += "," + col[1]

  
for k in sorted(csr, key = lambda x: float(x)):
    out.write(k + csr[k] + "\n")
out.close()
