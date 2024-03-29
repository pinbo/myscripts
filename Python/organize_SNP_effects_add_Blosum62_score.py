#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Junli Zhang <zhjl86@gmail.com>
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
# Before running this script
## extract effects and subsets from a vcf file
# cat Chen-confidence-region-subset.vcf | java -jar ../SnpSift.jar filter "( QUAL > 100 )" | ../scripts/vcfEffOnePerLine.pl | java -jar ../SnpSift.jar extractFields - CHROM POS REF ALT QUAL DP MQ "EFFECT" "GENE" "FEATUREID" "IMPACT" "HGVS_P" "GEN[*]" > Chen-confidence-region-subset-Qual100-SNPeff.txt
## convert GT to SNP table
# convert_vcf_calls_to_SNP_and_add_Blosum62_score.py Chen-confidence-region-subset-Qual100-SNPeff.txt Chen-confidence-region-subset-Qual100-SNPeff-converted.txt

# example
# ./organize_SNP_effects_add_Blosum62_score.py input_file output_file

# 2021-11-27 modifed to orgnize the effects to be one line for each position

# print usage
usage = '''
Example: organize_SNP_effects_add_Blosum62_score.py input_file output_file

To orgnize the effects to be one line for each position.

For output from command (no need GEN[*]):
cat confidence-region-subset.vcf | /path/to/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /path/to/SnpSift.jar extractFields - CHROM POS REF ALT QUAL DP MQ "EFFECT" "GENE" "FEATUREID" "IMPACT" "HGVS_P" "GEN[*]" > confidence-region-subset-SNPeff.txt

'''
print(usage)
### Imported
import sys
import re

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
	"*" : "*"
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
	"*" : "*"
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
# this header is only for commands of snpEffect
# cat your.vcf | pathto/scripts/vcfEffOnePerLine.pl | java -jar pathto/SnpSift.jar extractFields - CHROM POS REF ALT QUAL DP MQ "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].IMPACT" "ANN[*].HGVS_P" > effect-only-your.txt
header = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "MQ", "EFFECT", "GENE", "FEATUREID", "IMPACT", "HGVS_P", "AA_change", "Blosum62"]
out.write("\t".join(header) + "\n")

def gt2snp(allele_list, gt):
	c = "" # snp calls
	if "./." in gt: # missing
		c = "N"
	else:
		a, b = gt.split("/")
		if a == b: # homozygous
			c = allele_list[int(a)]
		else:
			c = "H" # heterozygous
	return c

pos_dict = {} # store positions
last_pos = [] # store infor for last position
with open(in_file) as file_one:
	header_line = next(file_one)
	# out.write(header_line + "\n")
	for line in file_one:
		# log.write(line)
		line = line.strip()
		#print("line number: ", n);
		# calculate BLOSUM62 score
		if line:
			if "stream" in line: # skip downstream and upstream variants
				continue
			ll = line.split("\t")
			B62 = ""
			AAsingle = ""
			if ll[11]:
				AAchange = ll[11].lstrip("p.") # p.Arg25His
				if "del" in AAchange or "ins" in AAchange or "*" in AAchange or "fs" in AAchange:
					AAsingle = AAchange
					B62 = "-10"
				else:
					try:
						ref_AA, pos, alt_AA = re.split(r'(\d+)', AAchange)
					except ValueError:
							AAsingle = "NA"
							continue
					if ref_AA in AA3letter2 and alt_AA in AA3letter2:
						AAsingle = AA3letter2[ref_AA] + pos + AA3letter2[alt_AA]
					if ll[7] == "missense_variant":
						try:
							B62 = B62table[B62header2.index(ref_AA)][B62header2.index(alt_AA)]
						except ValueError:
							B62 = "NA"
			if ll[10] == "HIGH":
				B62 = "-10"
			if ll[1] in pos_dict: # already present
				last_pos[9]+= "," + ll[9]
				if ll[11]:
					last_pos[11]+= "," + ll[11]
					last_pos[12]+= "," + AAsingle
				if B62:
					last_pos[13]+= "," + str(B62)
			else:
				if last_pos:
					out.write("\t".join(last_pos) + "\n")
				pos_dict[ll[1]] = 0
				last_pos = ll[0:12] + [AAsingle, str(B62)]
			# convert GT to SNPs
			# GTs = [x.split(":")[0] for x in ll[12:]]
			# ref = ll[2]
			# alt = ll[3].split(",") # there may be more than one alternative alleles
			# alleles = [ref] + alt # this way, 0 will be ref, and 1, 2 ... will be alternative allele
			# SNPs = [gt2snp(alleles, x) for x in GTs]
			# out.write("\t".join(ll[0:12] + [AAsingle, str(B62)] + SNPs) + "\n")

# write the last record after all lines
out.write("\t".join(last_pos) + "\n")
out.close()