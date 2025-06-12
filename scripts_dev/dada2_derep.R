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
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]"),
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

cat("DADA2 version:", as.character(packageVersion("dada2")), "\n")

fnFs <- sort(list.files(opt$input_dir, pattern = "_R1_trimmed.fastq.gz$", full.names = TRUE))
sample.names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))

filtFs <- file.path(opt$output_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output_dir, paste0(sample.names, "_R_filt.fastq.gz"))

for (i in seq_along(sample.names)) {
  sample <- sample.names[i]
  cat("Processing sample:", sample, "\n")

  derepF <- dada2::derepFastq(filtFs[i])
  derepR <- dada2::derepFastq(filtRs[i])

  saveRDS(derepF, file = file.path(opt$output_dir, "derep_objects", paste0(sample, "_F_derep.rds")))
  saveRDS(derepR, file = file.path(opt$output_dir, "derep_objects", paste0(sample, "_R_derep.rds")))

  rm(derepF, derepR)
  gc()
}

