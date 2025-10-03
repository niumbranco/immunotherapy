#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

# Load ASV inference results
dadaFs_paths <- snakemake@input[["dadaFs"]]
dadaRs_paths <- snakemake@input[["dadaRs"]]

dadaFs <- lapply(dadaFs_paths, readRDS)
dadaRs <- lapply(dadaRs_paths, readRDS)

seqtab <- makeSequenceTable(dadaFs)

saveRDS(seqtab, snakemake@output[["seqtab"]])