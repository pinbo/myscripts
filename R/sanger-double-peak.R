#!/usr/bin/Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Please give at least one ab1 file name (input file).\n")
}
library(sangerseqR)

for (i in args){
	outfile = paste0("chromatogram-", sub("ab1", "pdf", i))

	# which(strsplit(a,"")[[1]]!=strsplit(b,"")[[1]])

	## Test my samples
	hetsangerseq <- readsangerseq(i)
	# Making Basecalls
	hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
	chromatogram(hetcalls, width = 80, height = 1.5, trim5 = 20, trim3 = 200, showcalls = "both", filename = outfile)

	# Pairwise alignment
	# pa <- pairwiseAlignment(hetcalls@primarySeq, hetcalls@secondarySeq,
	#                         type = "global-local")
	# writePairwiseAlignments(pa)
}
