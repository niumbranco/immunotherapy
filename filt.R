#!/usr/bin/env Rscript

#load packages 
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
  library(tidyverse)
})

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing trimmed FASTQ files"),
  make_option(c("-o", "--output_dir"), type="character", default="dada2_out", help="Output directory [default %default]"),
  make_option(c("--trunc_len_f"), type="integer", default=240, help="Forward read truncation length [default %default]"),
  make_option(c("--trunc_len_r"), type="integer", default=200, help="Reverse read truncation length [default %default]"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

fnFs <- sort(list.files(opt$input_dir, pattern="_1.clean.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(opt$input_dir, pattern="_2.clean.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- paste0(sample.names, "_16S")

filtFs <- file.path(opt$output_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output_dir, paste0(sample.names, "_R_filt.fastq.gz"))

cat("Reading files from:", opt$input_dir, "\n")
print(list.files(opt$input_dir))

filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                          truncLen = c(opt$trunc_len_f, opt$trunc_len_r),
                          maxEE = c(2,2),
                          truncQ=2,
                          rm.phix = TRUE,
                          compress = TRUE,
                          multithread = opt$threads)

write.table(filt_out, file = file.path(opt$output_dir, "filtering_summary.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
