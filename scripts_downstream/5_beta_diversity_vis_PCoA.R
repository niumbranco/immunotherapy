library(tidyverse)

# working directory # where to save the results
work_dir <- "/Users/smartinezarbas/Repositories/nium/friesland_campina/data"
setwd(work_dir)

metadata <- read_tsv(file.path(work_dir,"tables", "metadata.tsv"), show_col_types = FALSE)

# read tb
bray_tb <- read.table(file.path(work_dir, "tables", paste0("bray_dist_matrix.tsv")),
                        header = TRUE, row.names = 1, sep = "\t")
bray_dist <-  as.dist(as.matrix(bray_tb))

# Run PCoA on Bray-Curtis
pcoa_bray <- pcoa(bray_dist)

# Create data frame for plotting
pcoa_df <- pcoa_bray$vectors[, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(.,metadata, by = c("sample_id" = "sample"))  # merge with sample metadata

pcoa_df <- pcoa_df %>%
  mutate(time = as.factor(time))

pcoa_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = condition, shape = time)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA - Bray-Curtis - beta diversity",
    x = paste0("PCoA 1 (", round(pcoa_bray$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa_bray$values$Relative_eig[2] * 100, 1), "%)")
  )

pdf(file = file.path(work_dir, "visuals" ,"2_beta_diversity_bray_PCoA.pdf"),
    width = 12/2.54,
    height = 12/2.54)
print(pcoa_plot)
dev.off()

# save table with coordinates, and name the columns with the explained percentage

library(glue)

# Create dynamic name
axis1_label <- glue("PCoA 1 ({round(pcoa_bray$values$Relative_eig[1] * 100, 1)}%)")
axis2_label <- glue("PCoA 2 ({round(pcoa_bray$values$Relative_eig[2] * 100, 1)}%)")

# Rename the column
coord_tb <- pcoa_df %>%
  select(sample_id, Axis.1, Axis.2) %>% 
  rename(!!axis1_label := Axis.1,
         !!axis2_label := Axis.2)
  
write_tsv(coord_tb, file.path(work_dir, "tables" ,"PCoA_coordinates.tsv"))

# # Use vegan::envfit() to fit feature vectors (species) onto the ordination. This estimates how each taxon’s abundance correlates with the PCoA axes
# library(vegan)
# library(ggrepel)
# 
# # Assume `matrix` is samples x species (same used in `vegdist`)
# # and `pcoa_bray$vectors` contains sample scores
# 
# fit <- envfit(pcoa_bray$vectors, matrix, permutations = 99)
# 
# species_vectors <- scores(fit, display = "vectors") %>%
#   as.data.frame() %>%
#   rownames_to_column("species")
# 
# pvals <- fit$vectors$pvals
# species_vectors$pval <- pvals
# 
# # Top N taxa by magnitude (vector length)
# top_species <- species_vectors %>%
#   #filter(pval < 0.1) %>% 
#   mutate(magnitude = sqrt(Axis.1^2 + Axis.2^2)) %>%
#   arrange(desc(magnitude)) %>% 
#   slice_head(n = 10)
# 
# top_opposite <- species_vectors %>%
#   mutate(magnitude = sqrt(Axis.1^2 + Axis.2^2)) %>%
#   filter(Axis.1 < 0 | Axis.2 < 0) %>%
#   arrange(desc(magnitude)) %>% 
#   slice_tail(n = 10)
# 
# top_species <- rbind(top_species, top_opposite)
# 
# # plot
# # Base PCoA plot
# colnames(pcoa_df) <- c("sample_id", "PCoA1", "PCoA2", "chamner", "time", "condition", "group") 
# 
# p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = condition)) +
#   geom_point(size = 3, alpha = 0.7) +
#   theme_minimal() +
#   labs(title = "PCoA with Top Taxa Arrows",
#        x = "PCoA 1",
#        y = "PCoA 2")
# 
# # Add arrows for top taxa
# pcoa_load_plot <- p + geom_segment(data = top_species,
#                                    aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
#                                    arrow = arrow(length = unit(0.2, "cm")),
#                                    color = "turquoise") +
#   geom_text_repel(data = top_species,
#                   aes(x = Axis.1, y = Axis.2, label = species),
#                   size = 3, 
#                   max.overlaps = Inf, 
#                   color = "turquoise",
#                   segment.color = "grey50")
# 
# pdf(file = file.path(work_dir, "beta_diversity_PCOA_loads.pdf"),
#     width = 12/2.54,
#     height = 12/2.54)
# print(pcoa_load_plot)
# dev.off()