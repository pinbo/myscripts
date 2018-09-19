#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Reverse complement short sequences line by line
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
# getclip | ./reverse_complement_short_sequence_line_by_line_from_stdin.py

import sys

# get reverse complement sequence
def ReverseComplement(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

# get stdin
# borrowed from https://stackoverflow.com/questions/4429966/how-to-make-a-python-script-pipeable-in-bash
if __name__ == "__main__":
	for line in sys.stdin:
		line = line.strip()
		#sys.stderr.write("DEBUG: got line: " + line)
		sys.stdout.write(ReverseComplement(line) + "\n")
