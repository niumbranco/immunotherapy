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
  make_option(c("--silva_train_set"), type="character", help="Path to SILVA training set fasta"),
  make_option(c("--silva_species_set"), type="character", help="Path to SILVA species fasta (optional)", default=NULL)
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))


# Read ASV table
# Read ASV table
seqtab.nochim <- readRDS(file.path(opt$output_dir, "asv_table.rds"))

# Split ASVs into chunks
asvs <- colnames(seqtab.nochim)
chunk_size <- 5000
chunks <- split(asvs, ceiling(seq_along(asvs) / chunk_size))

# Classify each chunk
taxa_list <- list()
for (i in seq_along(chunks)) {
  cat("Classifying chunk", i, "of", length(chunks), "\n")
  subset_tab <- seqtab.nochim[, chunks[[i]], drop = FALSE]
  taxa_chunk <- assignTaxonomy(subset_tab, opt$silva_train_set, multithread = opt$threads)
  taxa_list[[i]] <- taxa_chunk
  gc()
}

# Combine all taxonomy chunks
taxa <- do.call(rbind, taxa_list)

# Optional species assignment
if (!is.null(opt$silva_species_set)) {
  cat("Assigning species (in chunks)...\n")

  # Reuse the same ASV chunks
  taxa_species_list <- list()
  for (i in seq_along(chunks)) {
    cat("Species assignment for chunk", i, "of", length(chunks), "\n")
    # Subset taxonomy by ASVs
    taxa_chunk <- taxa[chunks[[i]], , drop = FALSE]
    taxa_species_chunk <- addSpecies(taxa_chunk, opt$silva_species_set, verbose = TRUE)
    taxa_species_list[[i]] <- taxa_species_chunk
    gc()
  }

  # Combine all species-assigned chunks
  taxa <- do.call(rbind, taxa_species_list)
}

# Save taxonomy table
saveRDS(taxa, file.path(opt$output_dir, "taxonomy.rds"))

