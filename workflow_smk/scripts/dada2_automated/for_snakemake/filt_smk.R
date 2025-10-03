#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(dada2)
  library(tidyverse)
})

# Get inputs from Snakemake
fnFs <- c(snakemake@input[["read_1"]])  #forward reads 
fnRs <- c(snakemake@input[["read_2"]])  #reverse reads 
sample_names <- snakemake@params[["sample_names"]]
output_dir <- snakemake@params[["output_dir"]]
trunc_len_f <- snakemake@params[["trunc_len_f"]]
trunc_len_r <- snakemake@params[["trunc_len_r"]]
filt_summary_path <- snakemake@output[["summary"]]

# Output files
filtFs <- c(snakemake@output[["filt1"]])
filtRs <- c(snakemake@output[["filt2"]])

# Create output directory if needed
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Apply filterAndTrim
filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                          truncLen = c(trunc_len_f, trunc_len_r),
                          maxEE = c(2,2),
                          truncQ = 2,
                          rm.phix = TRUE,
                          compress = TRUE,
                          multithread = TRUE)

# Save summary
write.table(filt_out, file = filt_summary_path,
            sep = "\t", quote = FALSE, col.names = NA)