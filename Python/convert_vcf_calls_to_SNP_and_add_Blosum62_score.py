#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
#

# example
# ./add_Blosum62_score.py input_file output_file

### Imported
import sys

AA3letter = {
	"I" : "Ile",
	"L" : "Leu",
	"V" : "Val",
	"F" : "Phe",
	"M" : "Met",
	"C" : "Cys",
	"A" : "Ala",
	"G" : "Gly",
	"P" : "Pro",
	"T" : "Thr",
	"S" : "Ser",
	"Y" : "Tyr",
	"W" : "Trp",
	"Q" : "Gln",
	"N" : "Asn",
	"H" : "His",
	"E" : "Glu",
	"D" : "Asp",
	"K" : "Lys",
	"R" : "Arg",
	"*" : "Stop"
}
AA3letter2 = {
	"Ile" : "I",
	"Leu" : "L",
	"Val" : "V",
	"Phe" : "F",
	"Met" : "M",
	"Cys" : "C",
	"Ala" : "A",
	"Gly" : "G",
	"Pro" : "P",
	"Thr" : "T",
	"Ser" : "S",
	"Tyr" : "Y",
	"Trp" : "W",
	"Gln" : "Q",
	"Asn" : "N",
	"His" : "H",
	"Glu" : "E",
	"Asp" : "D",
	"Lys" : "K",
	"Arg" : "R",
	"Stop" : "*"
}
# BLOSUM62
B62header = ["C", "S", "T", "P", "A", "G", "N", "D", "E", "Q", "H", "R", "K", "M", "I", "L", "V", "F", "Y", "W", "*"]
B62header2 = ["Cys", "Ser", "Thr", "Pro", "Ala", "Gly", "Asn", "Asp", "Glu", "Gln", "His", "Arg", "Lys", "Met", "Ile", "Leu", "Val", "Phe", "Tyr", "Trp", "Stop"]

B62table = [[9, -1, -1, -3, 0, -3, -3, -3, -4, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2, -4],
[-1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, -1, 0, -1, -2, -2, -2, -2, -2, -3, -4],
[-1, 1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, 0, -1, -2, -2, -2, -2, -2, -3, -4],
[-3, -1, 1, 7, -1, -2, -1, -1, -1, -1, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4, -4],
[0, 1, -1, -1, 4, 0, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -2, -2, -2, -3, -4],
[-3, 0, 1, -2, 0, 6, -2, -1, -2, -2, -2, -2, -2, -3, -4, -4, 0, -3, -3, -2, -4],
[-3, 1, 0, -2, -2, 0, 6, 1, 0, 0, -1, 0, 0, -2, -3, -3, -3, -3, -2, -4, -4],
[-3, 0, 1, -1, -2, -1, 1, 6, 2, 0, -1, -2, -1, -3, -3, -4, -3, -3, -3, -4, -4],
[-4, 0, 0, -1, -1, -2, 0, 2, 5, 2, 0, 0, 1, -2, -3, -3, -3, -3, -2, -3, -4],
[-3, 0, 0, -1, -1, -2, 0, 0, 2, 5, 0, 1, 1, 0, -3, -2, -2, -3, -1, -2, -4],
[-3, -1, 0, -2, -2, -2, 1, 1, 0, 0, 8, 0, -1, -2, -3, -3, -2, -1, 2, -2, -4],
[-3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3, -4],
[-3, 0, 0, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5, -1, -3, -2, -3, -3, -2, -3, -4],
[-1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5, 1, 2, -2, 0, -1, -1, -4],
[-1, -2, -2, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 1, 0, -1, -3, -4],
[-1, -2, -2, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4, 3, 0, -1, -2, -4],
[-1, -2, -2, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4, -1, -1, -3, -4],
[-2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6, 3, 1, -4],
[-2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7, 2, -4],
[-2, -3, -3, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2, 11, -4],
[-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]

### read files
in_file = sys.argv[1]
out_file = sys.argv[2]
out = open(out_file, "w")
log = open("log.txt", "w")
header = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "MQ", "ANN[*].EFFECT", "ANN[*].GENE", "ANN[*].FEATUREID", "ANN[*].IMPACT", "ANN[*].HGVS_P", "Blosum62", "2045A", "AvocetS", "Bakahtawar94", "Berkut", "CAP2", "D447-DW1", "DD", "DV92", "Excalibur", "G25-Yr15", "G3116", "Gredho", "Hahn-RW", "Inayama", "Kronos0", "Kronos280-1-Yr15", "LR23", "LR3", "Langdon", "PBW343", "PI347731", "PI503555", "PI519805", "Pavon", "RAC875", "RIL143", "RSI5", "Rusty", "Svevo", "UC1036", "UC1419", "Zavitan"]

out.write("\t".join(header) + "\n")

def gt2snp(allele_list, gt):
	c = "" # snp calls
	if "." in gt: # missing "./." or "."
		c = "N"
	else:
		a, b = gt.split("/")
		if a == b: # homozygous
			c = allele_list[int(a)]
		else:
			c = "H" # heterozygous
	return c

with open(in_file) as file_one:
	header_line = next(file_one)
	for line in file_one:
		log.write(line)
		line = line.strip()
		# calculate BLOSUM62 score
		if line:
			ll = line.split("\t")
			if ll[7] == "missense_variant":
				AAchange = ll[11] # p.Arg25His
				ref_AA = AAchange[2:5]
				alt_AA = AAchange[-3:]
				#pos = AAchange[5:-3]
				B62 = B62table[B62header2.index(ref_AA)][B62header2.index(alt_AA)]
				#out.write(line + "\t" + str(B62) + "\n")
			else:
				#out.write(line + "\n")
				B62 = "NA"
		# convert GT to SNPs
		GTs = ll[12:] # all GT calls
		ref = ll[2]
		alt = ll[3].split(",") # there may be more than one alternative alleles
		alleles = [ref] + alt # this way, 0 will be ref, and 1, 2 ... will be alternative allele
		SNPs = [gt2snp(alleles, x) for x in GTs]
		out.write("\t".join(ll[0:12] + [str(B62)] + SNPs) + "\n")

out.close()
log.close()