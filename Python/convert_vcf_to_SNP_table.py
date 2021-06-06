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
# Before running this script
## extract effects and subsets from a vcf file
# cat Chen-confidence-region-subset.vcf | java -jar ../SnpSift.jar filter "( QUAL > 100 )" | ../scripts/vcfEffOnePerLine.pl | java -jar ../SnpSift.jar extractFields - CHROM POS REF ALT QUAL DP MQ "EFFECT" "GENE" "FEATUREID" "IMPACT" "HGVS_P" "GEN[*]" > Chen-confidence-region-subset-Qual100-SNPeff.txt
## convert GT to SNP table
# convert_vcf_calls_to_SNP_and_add_Blosum62_score.py Chen-confidence-region-subset-Qual100-SNPeff.txt Chen-confidence-region-subset-Qual100-SNPeff-converted.txt

# example
# ./convert_vcf_calls_to_SNP_and_add_Blosum62_score.py input_file output_file

### Imported
import sys
import re

### read files
in_file = sys.argv[1]
out_file = sys.argv[2]
out = open(out_file, "w")
# ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DB_001", "DB_002", "DB_003", "DB_004", "DB_005", "DB_006", "DB_007", "DB_008", "DB_009", "DB_010", "DB_011", "DB_012", "DB_013", "DB_014", "DB_015DB_016", "DB_017", "DB_018", "DB_019", "DB_020", "DB_021", "DB_022", "DB_023", "DB_024", "DB_025", "DB_026", "DB_027", "DB_028", "DB_029", "DB_030", "Undetermined"]]

# out.write("\t".join(header) + "\n")

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
		line = line.strip()
		if line:
			ll = line.split("\t")
			if line.startswith("#"):
				if line.startswith("#CHROM"):
					out.write("\t".join(ll[0:2] + ll[3:6] + ["AC", "AN", "DP", "MQ"] + ll[9:]) + "\n")
					continue
				else:
					continue
			# get AC, AN, and DP and MQ
			info = dict(x.split("=") for x in ll[7].split(";"))
			AC = info["AC"]
			AN = info["AN"]
			DP = info["DP"]
			MQ = info["MQ"]
			# convert GT from number to alleles
			GTs = [x.split(":")[0] for x in ll[9:]]
			ref = ll[3]
			alt = ll[4].split(",") # there may be more than one alternative alleles
			alleles = [ref] + alt # this way, 0 will be ref, and 1, 2 ... will be alternative allele
			SNPs = [gt2snp(alleles, x) for x in GTs]
			# write output
			out.write("\t".join(ll[0:2] + ll[3:6] + [AC, AN, DP, MQ] + SNPs) + "\n")
out.close()
