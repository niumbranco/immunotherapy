#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(tidyverse)
})

# Define options
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing trimmed FASTQ files"),
  make_option(c("-o", "--output_dir"), type="character", default="dada2_out", help="Output directory [default %default]"),
  make_option(c("--trunc_len_f"), type="integer", default=240, help="Forward read truncation length [default %default]"),
  make_option(c("--trunc_len_r"), type="integer", default=200, help="Reverse read truncation length [default %default]"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]"),
  make_option(c("--silva_train_set"), type="character", help="Path to SILVA training set fasta"),
  make_option(c("--silva_species_set"), type="character", help="Path to SILVA species fasta (optional)", default=NULL)
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

# Define filtered output paths
fnFs <- sort(list.files(opt$input_dir, pattern = "_R1_trimmed.fastq.gz$", full.names = TRUE))
sample.names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))
filtFs <- file.path(opt$output_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# Error rates
errF <- learnErrors(filtFs, multithread = opt$threads)
errR <- learnErrors(filtRs, multithread = opt$threads)

saveRDS(errF, file.path(opt$output_dir, "errF.rds"))
saveRDS(errR, file.path(opt$output_dir, "errR.rds"))
