#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(dada2)
})

# Get inputs and outputs from Snakemake
filtFs <- snakemake@input[["filt_F"]]
filtRs <- snakemake@input[["filt_R"]]
errF_path <- snakemake@output[["errF"]]
errR_path <- snakemake@output[["errR"]]

cat("Learning error rates from all filtered FASTQ files...\n")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

saveRDS(errF, errF_path)
saveRDS(errR, errR_path)