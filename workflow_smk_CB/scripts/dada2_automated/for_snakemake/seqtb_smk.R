#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

# --- Load input paths ---
dadaFs_paths <- snakemake@input[["dadaFs"]]
dadaRs_paths <- snakemake@input[["dadaRs"]]

# --- Extract sample names from filenames ---
# Example: sample_1_dadaFs.rds → sample_1
get_sample_name <- function(path) {
  sub("_dadaFs\\.rds$", "", basename(path))
}

samples <- sapply(dadaFs_paths, get_sample_name)

# --- Load dada objects ---
dadaFs <- lapply(dadaFs_paths, readRDS)
dadaRs <- lapply(dadaRs_paths, readRDS)

# --- Assign names to each list element ---
names(dadaFs) <- samples
names(dadaRs) <- samples

# --- Build sequence table ---
seqtab <- makeSequenceTable(dadaFs)

# --- Assign sample names to rows ---
rownames(seqtab) <- samples

# --- Save ---
saveRDS(seqtab, snakemake@output[["seqtab"]])
