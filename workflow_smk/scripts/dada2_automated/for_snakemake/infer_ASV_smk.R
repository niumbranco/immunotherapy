#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(dada2)
})

# Load filtered FASTQ file paths (same logic as in err.R)
filtF <- snakemake@input[["filtF"]]
filtR <- snakemake@input[["filtR"]]
errF <- readRDS(snakemake@input[["errF"]])
errR <- readRDS(snakemake@input[["errR"]])

# Get outputs from Snakemake
dadaFs_out <- snakemake@output[["dadaF"]]
dadaRs_out <- snakemake@output[["dadaR"]]
mergers_out <- snakemake@output[["mergers"]]

# Dereplication
derepF <- derepFastq(filtF)
derepR <- derepFastq(filtR)

# inference
dadaFs <- dada(derepF, err=errF, multithread = TRUE)
dadaRs <- dada(derepR, err=errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepF, dadaRs, derepR)

saveRDS(dadaFs, dadaFs_out)
saveRDS(dadaRs, dadaRs_out)
saveRDS(mergers, mergers_out)
