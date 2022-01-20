#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  blast pangenome and extrat sequences
#  Copyright 2022 Junli Zhang <zhjl86@gmail.com>
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

#########################
def main(args):
	usage = '''
	This tool blasts the pangenome and extract sequences for later alignments.
	Please provide 4 arguments:
	  1. query fasta file (a few hundred bp that starts from ATG and is enough to differ homeologs);
	  2. output file name (no spaces allowed);
	  3. 5' extend length: for extracting promoter region of a gene. Usually 2kb is enough;
	  4. 3' extend length: usually your gene length.
	'''
	print (usage)
	query = args[1] # a  fasta file for the beginning of genes
	outfile = args[2]
	fivePrimeExtend = int(args[3])
	threePrimeExtend = int(args[4])
	reference = {
		"Svevo" : "/var/www/html/lab_blast/db/nucleotide/160802_Svevo_v2_pseudomolecules.fasta",
		"Julius" : "/var/www/html/lab_blast/db/nucleotide/170807_julius_MAGIC3_pseudomolecules.fasta",
		"CDC_Landmark" : "/var/www/html/lab_blast/db/nucleotide/170831_Landmark_pseudomolecules.fasta",
		"Jagger" : "/var/www/html/lab_blast/db/nucleotide/180529_Jagger_pseudomolecule_v1.1.fasta",
		"ArinaLrFor" : "/var/www/html/lab_blast/db/nucleotide/180808_ArinaLrFor_pseudomolecules_v3.fasta",
		"CDC_Stanley" : "/var/www/html/lab_blast/db/nucleotide/180902_Stanley_pseudomolecules_v1.2.fasta",
		"SY_Mattis" : "/var/www/html/lab_blast/db/nucleotide/181016_SY_Mattis_pseudomolecule_v1.fasta",
		"Lancer" : "/var/www/html/lab_blast/db/nucleotide/181120_lancer_pseudomolecule_v1.0.fasta",
		"Mace" : "/var/www/html/lab_blast/db/nucleotide/181120_mace_pseudomolecule_v1.0.fasta",
		"Norin61" : "/var/www/html/lab_blast/db/nucleotide/190307_Norin61_pseudomolecule_v1.1.fasta",
		"PI190962-Spelt" : "/var/www/html/lab_blast/db/nucleotide/190524_spelt_pseudomolecules_v1.0.fasta",
		"Fielder" : "/var/www/html/lab_blast/db/nucleotide/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta",
		"Claire" : "/var/www/html/lab_blast/db/nucleotide/Triticum_aestivum_Claire_EIv1.1.fa",
		"Robigus" : "/var/www/html/lab_blast/db/nucleotide/Triticum_aestivum_Robigus_EIv1.1.fa",
		"Cadenza" : "/var/www/html/lab_blast/db/nucleotide/Cadenza_EI_v1_arm-classified.fasta",
		"Kronos" : "/var/www/html/lab_blast/db/nucleotide/Kronos_EI_v1.fasta",
		"Paragon" : "/var/www/html/lab_blast/db/nucleotide/Paragon_EI_v1.fasta",
		"CS-Triticum_4.0" : "/var/www/html/lab_blast/db/nucleotide/GCA_002220415.3_Triticum_4.0_genomic.fna",
		"W7984-synthetic-wheat" : "/var/www/html/lab_blast/db/nucleotide/w7984.meraculous.scaffolds.Mar28.fa"
	}
	
	#step 1: blast
	out = open(outfile, "w")
	for ref in reference:
		print ("\n\n\n=========  Blast " + ref + "  =========\n")
		cmd2 = 'blastn -task blastn -db ' + reference[ref] + ' -query ' + query + ' -outfmt 6 -max_target_seqs 1 -word_size 30 -num_threads 2'
		print ("Step 2: Blast command:\n", cmd2)
		# call(cmd2, shell=True)
		p = subprocess.run(cmd2.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		print( 'exit status:', p.returncode )
		print( 'stdout:', p.stdout.decode() )
		print( 'stderr:', p.stderr.decode() )
		if p.stdout.decode().strip():
			col = p.stdout.decode().strip().split('\t')
			entryName = col[1]
			sStart = int(col[8])
			sEnd = int(col[9])
			extractStart = sStart - fivePrimeExtend
			extractEnd = sEnd + threePrimeExtend
			strand = "plus"
			if sStart > sEnd:
				extractEnd = sStart + fivePrimeExtend
				extractStart = sEnd - threePrimeExtend
				strand = "minus"
			out.write(">" + ref + " " + entryName + " " +  str(extractStart) + "-" + str(extractEnd) + " " + strand + "\n")
			cmd3 = 'blastdbcmd -db ' + reference[ref] + ' -entry ' + entryName + ' -strand ' + strand +  ' -range ' + str(extractStart) + "-" + str(extractEnd) + " -outfmt \"%s\""
			print ("Step 3: Blast command:\n", cmd3)
			p = subprocess.run(cmd3.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			print( 'exit status:', p.returncode )
			# print( 'stdout:', p.stdout.decode() )
			print( 'stderr:', p.stderr.decode() )
			out.write(p.stdout.decode().strip('"\n') + "\n")

	out.close()
	print ("\n\n\n Sequences are extracted successfully!")
	return 0

if __name__ == '__main__':
	import sys, subprocess
	# from subprocess import call
	sys.exit(main(sys.argv))
