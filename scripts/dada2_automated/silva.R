#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

# Argument parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Input directory (containing seqtab_nochim.rds)"),
  make_option(c("-o", "--output_dir"), type = "character", default = "dada2_out", help = "Directory to save taxonomy output"),
  make_option(c("-r", "--ref_db"), type = "character", help = "Path to SILVA reference database fasta file"),
  make_option(c("-s", "--species_ref"), type = "character", help = "Path to SILVA species assignment file"),
  make_option(c("-t", "--threads"), type = "integer", default = 4, help = "Number of threads")
)

# Parse args
opt <- parse_args(OptionParser(option_list = option_list))

# Load the non-chimeric sequence table
seqtab.nochim <- readRDS(file.path(opt$input_dir, "seqtab_nochim.rds"))
cat("Loaded sequence table with dimensions:", dim(seqtab.nochim), "\n")

# Assign taxonomy
cat("Assigning taxonomy using training set:", opt$ref_db, "\n")
taxa <- assignTaxonomy(seqtab.nochim, opt$ref_db, multithread = opt$threads)

# Assign species
cat("Assigning species using species assignment file:", opt$species_ref, "\n")
taxa <- addSpecies(taxa, opt$species_ref)

# Save results
output_file <- file.path(opt$output_dir, "taxa_silva.rds")
saveRDS(taxa, output_file)
cat("Taxonomy assignment saved to taxa_silva.rds\n")