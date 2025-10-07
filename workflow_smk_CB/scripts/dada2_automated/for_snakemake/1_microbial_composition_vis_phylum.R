
# upload packages
library(tidyverse)

# working directory # where to save the results
work_dir <- "/Users/smartinezarbas/Repositories/nium/friesland_campina/data"
setwd(work_dir)

relabs_in <- file.path(work_dir, paste0("tables/Relabs_P.tsv"))
relabs_tb <- read_tsv(relabs_in, show_col_types = F) %>% # spread table: rows name, columns sample
              gather(sample, value, -name)

# Add to "others" taxa catagory those elements with less than 0.5% relative abundance
relabs_tb <- relabs_tb %>%
  mutate(name = ifelse(value < 0.5, "others", name))

taxa_list <- relabs_tb$name %>% unique()

relabs_tb <- relabs_tb %>% 
  # factor taxa values, to keep the order (alphabetic and with "others" as last)
  mutate(name = factor(name, levels = c(setdiff(taxa_list, "others"), "others"))) 

# Generate a custom palette for the current phylum, excluding "others"
set.seed(12)  # Ensures reproducible colors
color_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(taxa_list) - 1)

# Append "darkgrey" for "others"
custom_colors <- c(color_palette, "darkgrey")

legend_labels <- setNames(taxa_list, names(custom_colors))
legend_labels["others"] <- "others: group of < 0.5% taxa"

## split samples ID, as it has metadata information: chamber, time point, condition
relabs_tb <- relabs_tb %>% ungroup() %>% 
  mutate(
    chamber = str_extract(sample, "(?<=C)\\d"),
    time = str_extract(sample, "(?<=T)\\d+"),
    condition = str_extract(sample, "[FMH]")  # assuming always one of these
  ) %>% 
  mutate(group = paste0("t", time, "_", condition)) %>% 
  mutate(
    time = as.character(time),  # ensure it's not numeric
    time = str_trim(time),      # remove whitespace
    time = fct_relevel(time, "0", "6", "12", "24")  # reorder as factor
  )

# now, plot by groups
conditions <- c("M", "H", "F")

plots_list <- function(relabs_tb, cond_i){
  plot <- relabs_tb %>% 
    mutate(time = factor(as.character(time), levels = c("0", "6", "12", "24"))) %>%  # relevel here
    filter(condition == cond_i) %>% 
    ggplot(aes(x = sample, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors, labels = legend_labels) +
    labs(
      title = paste0("Condition: ", cond_i),
      x = "Sample",
      y = "Abundance",
      fill = "Phylum"
    ) +
    coord_cartesian(ylim = c(0, 100)) +  # better than `ylim()` to avoid dropped bars
    facet_wrap(~ time, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) 
  
  return(plot)
}

plot_M <- plots_list(relabs_tb, "M") +
  ylab("Condition: M") +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank()) +
  labs(title = "Microbial composition: phylum")

plot_F <- plots_list(relabs_tb, "F") + 
  labs(title = "") +
  ylab("Condition: F") +
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_blank())

plot_H <- plots_list(relabs_tb, "H") +
  labs(title = "") +
  ylab("Condition: H") 

#library(patchwork)
comb_plot <- plot_M + plot_F +  plot_H + plot_layout(ncol = 1)
print(comb_plot)

pdf(file = file.path(work_dir, "visuals", "0_microbial_composition_conditions_phylum.pdf"),
    height = 20/2.54,
    width = 28/2.54)
print(comb_plot)
dev.off()
