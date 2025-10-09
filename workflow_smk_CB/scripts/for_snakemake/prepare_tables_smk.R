#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# --- Inputs and outputs from Snakemake ---
seqtab_file <- snakemake@input[["seqtab"]]
taxa_file <- snakemake@input[["taxa"]]
out_phylum <- snakemake@output[["phylum"]]
out_genus  <- snakemake@output[["genus"]]

# --- Load DADA2 objects ---
seqtab <- readRDS(seqtab_file)
taxa <- readRDS(taxa_file)

# Ensure the ASV names match
common_asv <- intersect(colnames(seqtab), rownames(taxa))
seqtab <- seqtab[, common_asv]
taxa <- taxa[common_asv, ]

cat(sprintf("Keeping %d ASVs common to both seqtab (%d) and taxa (%d)\n",
            length(common_asv), ncol(seqtab), nrow(taxa)))

# Convert ASV abundance table to long format
seqtab_df <- as.data.frame(t(seqtab)) %>%
  tibble::rownames_to_column("ASV")

# Combine taxonomy info with counts
taxa_df <- as.data.frame(taxa, check.names = FALSE)
colnames(taxa_df) <- make.unique(colnames(taxa_df))
taxa_df <- tibble::rownames_to_column(taxa_df, "ASV")

merged_df <- left_join(taxa_df, seqtab_df, by = "ASV")

# Define ranks to summarize
ranks <- list(
  Phylum = out_phylum,
  Genus = out_genus
)

# --- Summarize counts per rank ---
for (r in names(ranks)) {
  message("Summarizing counts at rank: ", r)
  
  # Skip if column missing (e.g. not all ASVs classified to Species)
  if (!r %in% colnames(merged_df)) {
    warning(paste("Column", r, "not found in taxonomy table; skipping."))
    next
  }
  
  rank_counts <- merged_df %>%
    filter(!is.na(.data[[r]])) %>%
    group_by(.data[[r]]) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    rename(name = !!r)
  
  # Write to TSV
  write_tsv(rank_counts, ranks[[r]])

# --- Relative Abundances ---
message("  • Calculating relative abundances")

rank_relab <- rank_counts %>%
  mutate(across(where(is.numeric), \(x) ifelse(is.na(x), 0, x))) %>%
  gather(sample, counts, -name) %>%
  group_by(sample) %>%
  mutate(relab = counts / sum(counts)) %>%
  select(-counts) %>%
  spread(sample, relab)

relab_file <- sub("taxonomy_counts_", "Relabs_", ranks[[r]])
write_tsv(rank_relab, relab_file)

# --- CLR transformation ---
message("  • Calculating CLR transformation")

# Convert to matrix and replace zeros with small pseudocount
rank_matrix <- rank_counts %>%
  column_to_rownames("name") %>%
  as.matrix()
rank_matrix[rank_matrix == 0] <- 1e-6

clr_transf <- compositions::clr(rank_matrix) %>%
  as_tibble(rownames = "name")

clr_file <- sub("taxonomy_counts_", "CLR_", ranks[[r]])
write_tsv(clr_transf, clr_file)
}

message("\n✅ Taxonomic, relative abundance, and CLR tables successfully written.")
