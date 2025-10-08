#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(dada2)
})

# Snakemake I/O
seqtab_file <- snakemake@input[["seqtab"]]
out_file <- snakemake@output[["seqtab_nochim"]]

# Load the sequence table
seqtab <- readRDS(seqtab_file)

# Remove chimeras using consensus method (recommended by DADA2 authors)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

# Preserve rownames
rownames(seqtab_nochim) <- rownames(seqtab)

# Save the chimera-free table
saveRDS(seqtab_nochim, out_file)

# Optional summary message
cat(sprintf("Chimera removal complete. %d sequences remaining out of %d (%.2f%% retained)\n",
            ncol(seqtab_nochim), ncol(seqtab),
            (ncol(seqtab_nochim) / ncol(seqtab)) * 100))