#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  convert "primer-name sequence" to fasta format for primer blast. Get input from stdin
#  
#  Copyright 2018 Junli Zhang <zhjl86@gmail.com>
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

# example
# getclip | extract-sequence-pairwise-blast.py

import sys

# get stdin
# borrowed from https://stackoverflow.com/questions/4429966/how-to-make-a-python-script-pipeable-in-bash
if __name__ == "__main__":
	n = 0
	outseq = ">extracted_seq\n"
	for line in sys.stdin:
		n += 1 # line number
		line = line.strip()
		if line: # if not blank
			#print line
			if n % 4 == 3: # subject sequenc
				seq = line.split()[2]
				outseq += seq
	#sys.stderr.write("DEBUG: got line: " + line)
	sys.stdout.write(outseq + "\n")
