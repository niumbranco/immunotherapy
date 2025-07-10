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
  make_option(c("-t", "--threads"), type="integer", default=2, help="Number of threads [default %default]"),
  make_option(c("-f", "--suffix_1"), type="character", default="_F_filt.fastq.gz", help="Suffix of clean-trimmered forward reads [default %default]"),
  make_option(c("-r", "--suffix_2"), type="character", default="_R_filt.fastq.gz", help="Suffix of clean-trimmered reverse reads [default %default]")
  )

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

# Create output dir
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# List input files
fnFs <- sort(list.files(opt$input_dir, pattern = paste0(opt$suffix_1, "$"), full.names = TRUE))
fnRs <- sort(list.files(opt$input_dir, pattern = paste0(opt$suffix_2, "$"), full.names = TRUE))
sample.names <- gsub(opt$suffix_1, "", basename(fnFs))

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
