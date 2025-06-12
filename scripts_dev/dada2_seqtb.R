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
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]")
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

fnFs <- sort(list.files(opt$input_dir, pattern = "_R1_trimmed.fastq.gz$", full.names = TRUE))
sample.names <- gsub("_R1_trimmed.fastq.gz", "", basename(fnFs))

filtFs <- file.path(opt$output_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(opt$output_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# error objects
errF <- readRDS(file.path(opt$output_dir, "errF.rds"))
errR <- readRDS(file.path(opt$output_dir, "errR.rds"))

# Initialize empty lists to collect results
dadaFs <- list()
dadaRs <- list()
mergers <- list()

for (sample in sample.names) {
  cat("Inferring sample:", sample, "\n")

  # Load dereplicated objects
  derepF <- readRDS(file.path(opt$output_dir, "derep_objects", paste0(sample, "_F_derep.rds")))
  derepR <- readRDS(file.path(opt$output_dir, "derep_objects", paste0(sample, "_R_derep.rds")))

  # Sample inference
  dadaFs[[sample]] <- dada(derepF, err=errF, multithread=opt$threads)
  dadaRs[[sample]] <- dada(derepR, err=errR, multithread=opt$threads)

  # Merge paired reads
  mergers[[sample]] <- mergePairs(dadaFs[[sample]], derepF, dadaRs[[sample]], derepR)

  # Clean up memory
  rm(derepF, derepR)
  gc()
}

# Create sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = opt$threads)

# Save outputs
saveRDS(seqtab, file.path(opt$output_dir, "seqtab.rds"))
saveRDS(seqtab.nochim, file.path(opt$output_dir, "asv_table.rds"))

