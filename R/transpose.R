#!/usr/bin/Rscript --vanilla
# Function: transpose a dataset
# Time: 6/1/2015
# Author: Junli Zhang <zhjl86@gmail.com>
# example: transpose.R file1.txt
# input order is csr order (serpentine by column)
#cat("\nInput is a file for transposing !!!\n\n")
args <- commandArgs(trailingOnly = TRUE)
#nrow <- as.numeric(args[1])
filename <- args[1]
dd <- read.table(filename)
dd = t(dd)
write.table(dd, stdout(), sep="\t", row.names=F, col.names=F, quote=F)
