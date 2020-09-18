#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Junli Zhang <zhjl86@gmail.com>
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

## This script can be used to demultiplex samples user mixed as one sample
## these samples were preapred with two round PCRs
## 1st round:
# 5’-TCCTCTGTCACGGAAGCG-Your-forward-primer -3’
# 5’-TTTAGCCTCCCCACCGAC-Your-reverse-primer -3
## 2nd round: add barcodes (XXXXXXXX and YYYYYYYY are 8bp barcodes)
# 5’-XXXXXXXXTCCTCTGTCACGGAAGCG -3’
# 5’-YYYYYYYYTTTAGCCTCCCCACCGAC -3
## Then mix them and submit as one sample for NGS sequencing
## The sequencing facility then prepare libraries by adding their sequencing adapters and barcodes
## the comeback results have only our PCR products: 
## XXXXXXXXTCCTCTGTCACGGAAGCG... to 150 bp
## YYYYYYYYTTTAGCCTCCCCACCGAC... to 150 bp

## This script separate samples by XXXXXXXX and YYYYYYYY combinations

import sys

# classes
class sample(object):
	def __init__(self):
		self.name = ''
		self.R1 = [] # 4 lines of read1
		self.R2 = [] # 4 lines of read2

# function to process barcode file : 3 columns: sample name, left barcode, right barcode
def get_barcode(infile):
	barcodes = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line:
				name, leftbarcode, rightbarcode = line.split()
				barcodes[leftbarcode[-5:] + "-" + rightbarcode[-5:]] = name
	return barcodes

## arguments
barcodeFile = sys.argv[1] # 3 columns: sample name, left barcode, right barcode
fastqFile = sys.argv[2] # fastq file: both R1 and R2 are in the same file next to each other

## demultiplex
dictBarcode = get_barcode(barcodeFile)
dictSample = {}
n1 = 0
n2 = 1
sample_name = ""
leftAdapter  = "TCCTCT" # random_adapter_F TCCTCTGTCACGGAAGCG
rightAdapter = "TTTAGC" # random_adapter_R TTTAGCCTCCCCACCGAC
# I found the sequencing output sometimes has incomplete barcode.
# The right 5 bps are already unique for the random 8-bp barcodes I designed
# so I only use the right 5 bp here
barcodeLenToCheck = 5 
with open(fastqFile) as file_one:
	for line in file_one:
		line = line.strip()
		if line.startswith("@") and " " in line: # I found sometimes quality line (the 4th line) also starts with @, but they have no space
		# example: @M02850:171:000000000-J54GV:1:1101:19874:1325 1:N:0:1
			#print(line + "\n")
			sample_name, read = line.split(" ") # read is R1 or R2, here 1:N:0:1
			n1 = 2 # counter of the 4 lines of each read
			if sample_name in dictSample: # should be R2
				dictSample[sample_name].R2 = [line]
				n2 = 2 # now reads are R2
			else:
				ss = sample()
				ss.R1 = [line]
				dictSample[sample_name] = ss
		elif n2 == 1: # input to R1
			dictSample[sample_name].R1.append(line)
			n1 += 1
		elif n2 == 2 and n1 < 4:
			dictSample[sample_name].R2.append(line)
			n1 += 1
		elif n2 == 2 and n1 == 4: # ready to write
			dictSample[sample_name].R2.append(line)
			n2 = 1 # reset to 1
			leftbarcode = ""
			rightbarcode = ""
			R1First20bp = dictSample[sample_name].R1[1][:20] # the first 20 bp of R1
			R2First20bp = dictSample[sample_name].R2[1][:20] # the first 20 bp of R2
			P1 = R1First20bp.find(leftAdapter) # position of the left adpator in R1
			P2 = R2First20bp.find(leftAdapter)
			P3 = R1First20bp.find(rightAdapter)
			P4 = R2First20bp.find(rightAdapter)
			if P1 >= barcodeLenToCheck: # if there are still at least 5 bp on the left
				leftbarcode = R1First20bp[(P1 - barcodeLenToCheck):P1]
			if P2 >= barcodeLenToCheck:
				leftbarcode = R2First20bp[(P2 - barcodeLenToCheck):P2]
			if P3 >= barcodeLenToCheck:
				rightbarcode = R1First20bp[(P3 - barcodeLenToCheck):P3]
			if P4 >= barcodeLenToCheck:
				rightbarcode = R2First20bp[(P4 - barcodeLenToCheck):P4]
			barcode = leftbarcode + "-" + rightbarcode
			if barcode in dictBarcode:
				out = dictBarcode[barcode] + ".fastq"
			else:
				out = "unassigned.fastq"
			outfile = open(out, 'a')
			outfile.write('\n'.join(dictSample[sample_name].R1 + dictSample[sample_name].R2) + "\n")
			outfile.close()
			del dictSample[sample_name] # delete this entry to save memory

## write all entries kept in dictSample into unassigned.fastq
out = "unassigned2.fastq"
outfile = open(out, 'a')
for sample_name in dictSample:
	outfile.write('\n'.join(dictSample[sample_name].R1 + dictSample[sample_name].R2) + "\n")

outfile.close()
