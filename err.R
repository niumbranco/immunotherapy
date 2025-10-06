#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

# Define options
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing trimmed FASTQ files"),
  make_option(c("-o", "--output_dir"), type="character", default="dada2_out", help="Output directory [default %default]"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]")
) 

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Get filtered FASTQ files
filtFs <- sort(list.files(opt$input_dir, pattern = "_F_filt.fastq.gz$", full.names = TRUE))
filtRs <- sort(list.files(opt$input_dir, pattern = "_R_filt.fastq.gz$", full.names = TRUE))

# Error rates
errF_path <- file.path(opt$output_dir, "errF.rds")
errR_path <- file.path(opt$output_dir, "errR.rds")

if (!file.exists(errF_path) || !file.exists(errR_path)) {
  cat("Learning error rates...\n")
  errF <- learnErrors(filtFs, multithread = opt$threads)
  errR <- learnErrors(filtRs, multithread = opt$threads)
  saveRDS(errF, errF_path)
  saveRDS(errR, errR_path)
} else {
  cat("Loading existing error models\n")
  errF <- readRDS(errF_path)
  errR <- readRDS(errR_path)
}

# Plot and save error plots
# pdf(file.path(opt$output_dir, "errF_plot.pdf"))
# plotErrors(errF, nominalQ=TRUE)
# dev.off()

# pdf(file.path(opt$output_dir, "errR_plot.pdf"))
# plotErrors(errR, nominalQ=TRUE)
# dev.off() 
