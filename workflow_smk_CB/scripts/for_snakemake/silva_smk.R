#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

seqtab <- readRDS(snakemake@input[["seqtab"]])
silva_train <- snakemake@params[["silva_train"]]
silva_species <- snakemake@params[["silva_species"]]

ASV_size <- 100 #number of ASV that are treated at the same time
n_col <- ncol(seqtab)
ASV_part <- split(1:n_col, ceiling(seq_along(1:n_col)/ASV_size))

taxa_list <- list() 
row_counts <- c()

#loop
for (i in seq_along(ASV_part)) {
  message(sprintf("Processing ASV part %d of %d...", i, length(ASV_part)))
  cols <- ASV_part[[i]]
  sub_seqtab <- seqtab[, cols, drop = FALSE]
  
  #assign taxonomy
  taxa_asv_part <- tryCatch({
    taxa_asv_part <- assignTaxonomy(sub_seqtab, silva_train, multithread = TRUE)
    taxa_asv_part <- addSpecies(taxa_asv_part, silva_species)
    row_counts <- c(row_counts, nrow(taxa_asv_part))
    taxa_asv_part
  }, error = function(e) {
    warning(sprintf("Error in part %d: %s", i, e$message))
    NULL
  })
  
  if (!is.null(taxa_asv_part)) {
    taxa_list[[length(taxa_list) + 1]] <- taxa_asv_part
  }
}

# ensure consistent row counts
if (length(taxa_list) == 0) {
  stop("No valid taxa parts were generated.")
}

#pad with NA if necessary
max_rows <- max(sapply(taxa_list, nrow))
taxa_list <- lapply(taxa_list, function(x) {
  if (nrow(x) < max_rows) {
    pad <- matrix(NA, nrow = max_rows - nrow(x), ncol = ncol(x))
    colnames(pad) <- colnames(x)
    rbind(x, pad)
  } else {
    x
  }
})

# Assignation
taxa <- do.call(cbind, taxa_list)

saveRDS(taxa, snakemake@output[["taxa"]])