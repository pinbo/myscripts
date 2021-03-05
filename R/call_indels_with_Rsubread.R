## 2021-03-04
# Copyright 2021 Junli Zhang <zhjl86@gmail.com>
# This program is free software but WITHOUT ANY WARRANTY

## Install package Rsubread first
# https://bioconductor.org/packages/release/bioc/html/Rsubread.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")

# Windows users need to download the Windows Binary zip file at the bottom of the page
# unzip it and copy the "Rsubread" folder to the Rlib folder(use .libPaths() to find the folder position)

###################################################

## Step 1: create function "call_indels" by running the code below step 2

## Step 2: modify the paramters for your case and run "call_indels" to make bams and call indels

call_indels(
  fastq_dir = "../demultiplexed-samples/",  # the folder with all your demultiplexed fastq files
  ref_fa = "../reference.fasta",     # the fasta file with all your gene sequences
  # most of the time, you do not need to change the parameters below 
  read1_suffix = "_R1_001.fastq.gz", # fastq file name suffix, for files like "xxxx_R1_001.fastq.gz"
  read2_suffix = "_R2_001.fastq.gz",
  outdir = "results",                # output result folder names, can be relative path like "../result"
  detectSV = T # detect structural variants: indels > 15bp, inversions etc
  )





## The main function
# Run the code below first before calling call_indels
call_indels = function(fastq_dir, ref_fa, read1_suffix="_R1_001.fastq.gz", read2_suffix="_R2_001.fastq.gz", outdir="results", detectSV = TRUE){
  require("Rsubread")
  
  # build index for alignment
  index_prefix="reference_index"
  buildindex(basename=index_prefix, reference=ref_fa, indexSplit=T, memory = 2000)
  
  # alignment
  files = list.files(fastq_dir)
  read1_files = grep(read1_suffix, files, value = T)
  read2_files = gsub(read1_suffix, read2_suffix, read1_files)
  out_files = gsub(read1_suffix, ".bam", read1_files)
  # create a folder for bam files
  if (!dir.exists(outdir))  dir.create(outdir, recursive=T)
  
  dum = lapply(1:length(read1_files), function(i) {
         align(index=index_prefix,
               readfile1=file.path(fastq_dir, read1_files[i]),
               readfile2=file.path(fastq_dir, read2_files[i]),   
               input_format='gzFASTQ', 
               output_file=file.path(outdir, out_files[i]),
               indels = 15,
               detectSV = detectSV,
               sortReadsByCoordinates = T,
               type = 1, # DNA
               readGroupID = as.character(i),
               readGroup = paste0("SM:", gsub(".bam","", out_files[i])),
               nthreads=1) # a bug: have to be 1 if need sorted bam
               })
  
  ## merge indel calls
  indel_files = list.files(outdir, "indel.vcf")
  indels <- do.call("rbind", lapply(indel_files, FUN = function(file) {
    tryCatch({
      tmp <- read.table(file.path(outdir, file), comment.char="#", sep="\t")
      sampleID = gsub(".bam.indel.vcf", "", file)
      cbind(Sample=sampleID, tmp)
    }, error=function(e) NULL)
  }))
  if (!is.null(indels)){ # if there are indels
    names(indels) = c("Sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    #gsub("INDEL;DP=", "", indels$INFO)
    counts = data.frame(do.call('rbind', strsplit(gsub("INDEL;DP=", "", indels$INFO),';SR=',fixed=TRUE)))
    counts = data.frame(sapply(counts, as.numeric))
    names(counts) = c("TotalCoverage", "indelCoverage")
    counts$indelPercent = counts$indelCoverage/counts$TotalCoverage*100
    indels = cbind(indels[-c(4,8,9)], counts)
    indels$indelSize = nchar(indels$ALT) - nchar(indels$REF) 
  }
  ## merge big indels (>15 bp) and potential structure variations
  big_indel_files = list.files(outdir, "breakpoints.vcf")
  big_indels <- do.call("rbind", lapply(big_indel_files, FUN = function(file) {
    tryCatch({
      tmp <- read.table(file.path(outdir, file), comment.char="#", sep="\t")
      sampleID = gsub(".bam.breakpoints.vcf", "", file)
      cbind(Sample=sampleID, tmp)
    }, error=function(e) NULL)
  }))
  if (!is.null(big_indels)){
    names(big_indels) = c("Sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    big_indels$indelCoverage = as.numeric(gsub("SVTYPE=BND;MATEID=bnd_.*;SR=", "", big_indels$INFO))
    # split big_indels into 2 files
    oddn = seq(1, nrow(big_indels),2)
    oddrows = big_indels[oddn, ]
    evenrows = big_indels[oddn+1, ]
    newdata = cbind(oddrows[1:3],evenrows[3],oddrows[10])
    names(newdata) = c("Sample", "CHROM", "Start", "End", "indelCoverage")
    newdata$indelSize = newdata$End - newdata$Start
  }
  ## write summary file
  cat("Making summary files ...\n\n")
  write.csv(indels,  "Summary_indels_less_than_16bp.csv", row.names = F)
  write.csv(newdata, "Summary_indels_more_than_15bp.csv", row.names = F)
  cat("Summary_indels_less_than_16bp.csv and Summary_indels_more_than_15bp.csv are created.\n\n")
  
  ## remove the indexed reference, they are quite big
  cat("Removing indexed references ...\n\n")
  index_files = list.files(".", "reference_index")
  for (f in index_files) file.remove(f)
}

