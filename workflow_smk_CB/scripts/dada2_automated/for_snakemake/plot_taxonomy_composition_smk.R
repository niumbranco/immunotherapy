suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})

# --- Snakemake inputs and outputs ---
inputs  <- snakemake@input
outputs <- snakemake@output

# --- Helper function to make plots ---
plot_taxonomy <- function(in_file, per_sample_pdf, overview_pdf, rank_label) {
  
  # Read table
  relabs_tb <- read_tsv(in_file, show_col_types = FALSE) %>%
    gather(sample, value, -name)
  
  # Group rare taxa (<0.5%) into "others"
  relabs_tb <- relabs_tb %>%
    mutate(name = ifelse(value < 0.5, "others", name))
  
  # Set order of taxa
  taxa_list <- unique(relabs_tb$name)
  relabs_tb <- relabs_tb %>%
    mutate(name = factor(name, levels = c(setdiff(taxa_list, "others"), "others")))
  
  # Define color palette
  set.seed(12)
  color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(taxa_list) - 1)
  custom_colors <- c(color_palette, "darkgrey")
  
  # --- (1) Per-sample plot ---
  p1 <- ggplot(relabs_tb, aes(x = sample, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Microbial composition per sample (", rank_label, " level)", sep = ""),
      x = "Sample",
      y = "Relative abundance (%)",
      fill = rank_label
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(file = per_sample_pdf, height = 18/2.54, width = 28/2.54)
  print(p1)
  dev.off()
  
  # --- (2) Overview plot (mean across samples) ---
  overview_tb <- relabs_tb %>%
    group_by(name) %>%
    summarise(mean_value = mean(value, na.rm = TRUE)) %>%
    arrange(desc(mean_value))
  
  p2 <- ggplot(overview_tb, aes(x = reorder(name, -mean_value), y = mean_value, fill = name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = custom_colors) +
    labs(
      title = paste("Overall microbial composition (", rank_label, " level)", sep = ""),
      x = "Taxon",
      y = "Mean relative abundance (%)",
      fill = rank_label
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(file = overview_pdf, height = 14/2.54, width = 20/2.54)
  print(p2)
  dev.off()
}

# --- Generate both plots for Phylum & Genus ---
plot_taxonomy(
  inputs[["phylum"]],
  outputs[["phylum_pdf"]],
  outputs[["phylum_overview"]],
  "Phylum"
)

plot_taxonomy(
  inputs[["genus"]],
  outputs[["genus_pdf"]],
  outputs[["genus_overview"]],
  "Genus"
)
