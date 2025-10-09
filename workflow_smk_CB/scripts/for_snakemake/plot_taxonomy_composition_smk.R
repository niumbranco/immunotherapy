#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})

# --- Snakemake inputs and outputs ---
inputs  <- snakemake@input
outputs <- lapply(snakemake@output, normalizePath, mustWork = FALSE)

# --- Helper function to compute true overall abundance from counts ---
compute_overall_abundance <- function(counts_file) {
  counts_tb <- read_tsv(counts_file, show_col_types = FALSE)
  
  counts_long <- counts_tb %>%
    gather(sample, counts, -name) %>%
    group_by(name) %>%
    summarise(total_counts = sum(counts, na.rm = TRUE)) %>%
    mutate(global_relab = total_counts / sum(total_counts)) %>%
    arrange(desc(global_relab))
  
  return(counts_long)
}

# --- Main plotting function ---
plot_taxonomy <- function(relabs_file, counts_file, per_sample_pdf, overview_pdf, rank_label) {
  
  # (1) Load relative abundance table (0–1 scale)
  relabs_tb <- read_tsv(relabs_file, show_col_types = FALSE) %>%
    gather(sample, value, -name)
  
  # Group rare taxa (<0.005 = 0.5%) into "others"
  relabs_tb <- relabs_tb %>%
    mutate(name = ifelse(value < 0.005, "others", name))
  
  # Define consistent ordering and colors
  taxa_list <- unique(relabs_tb$name)
  relabs_tb <- relabs_tb %>%
    mutate(name = factor(name, levels = c(setdiff(taxa_list, "others"), "others")))
  
  set.seed(12)
  color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(taxa_list) - 1)
  custom_colors <- c(color_palette, "darkgrey")
  
  # --- (A) Per-sample stacked barplot (values 0–1) ---
  p1 <- ggplot(relabs_tb, aes(x = sample, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Microbial composition per sample (", rank_label, " level)", sep = ""),
      x = "Sample",
      y = "Relative abundance (0-1)",
      fill = rank_label
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    )
  
  # --- Wider output for better readability ---
  pdf(file = per_sample_pdf, height = 18/2.54, width = 38/2.54)
  print(p1)
  dev.off()
  
  # --- (B) Overall composition pie plot ---
  overview_tb <- compute_overall_abundance(counts_file)
  
  p2 <- ggplot(overview_tb, aes(x = "", y = global_relab, fill = name)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Overall microbial composition (", rank_label, " level)", sep = ""),
      fill = rank_label
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  pdf(file = overview_pdf, height = 16/2.54, width = 20/2.54)
  print(p2)
  dev.off()
}

# --- Generate both plots for Phylum & Genus ---
plot_taxonomy(
  relabs_file   = inputs[["phylum"]],
  counts_file   = inputs[["phylum_counts"]],
  per_sample_pdf = outputs[["phylum_pdf"]],
  overview_pdf   = outputs[["phylum_overview"]],
  rank_label = "Phylum"
)

plot_taxonomy(
  relabs_file   = inputs[["genus"]],
  counts_file   = inputs[["genus_counts"]],
  per_sample_pdf = outputs[["genus_pdf"]],
  overview_pdf   = outputs[["genus_overview"]],
  rank_label = "Genus"
)
