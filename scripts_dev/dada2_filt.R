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
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]")
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

# Create output dir
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# List input files
fnFs <- sort(list.files(opt$input_dir, pattern = "_R1_trimmed.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(opt$input_dir, pattern = "_R2_trimmed.fastq.gz$", full.names = TRUE))
sample.names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))

# Define filtered output paths
filtFs <- file.path(opt$output_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# Filtering
filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                          truncLen = c(opt$trunc_len_f, opt$trunc_len_r),
                          maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                          compress = TRUE, multithread = opt$threads)
# Save processing summary
write.table(filt_out, file.path(opt$output_dir, "filter_stats.tsv"), sep="\t", quote=FALSE)
