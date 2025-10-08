library(tidyverse)
library(vegan)
library(ape)

# working directory # where to save the results
work_dir <- "/Users/smartinezarbas/Repositories/nium/friesland_campina/data"
setwd(work_dir)

relabs_in <- file.path(work_dir, "tables", paste0("Relabs_S.tsv"))
relabs_tb <- read_tsv(relabs_in, show_col_types = F) %>% # spread table: rows name, columns sample
  gather(sample, value, -name) %>% 
  group_by(sample) %>% 
  mutate(relab = value*100/sum(value)) %>% 
  select(-value) %>% 
  mutate(value = relab) %>% select(-relab)
  
# samples are rows and taxa are columns
matrix <- relabs_tb %>% 
  spread(name, value) %>% 
  column_to_rownames("sample") %>% 
  as.matrix()

# save matrix
write.table(matrix, file.path(work_dir, "tables", paste0("Relabs_S_matrix.tsv")))

alpha_diversity <- data.frame(
  Sample_ID = rownames(matrix),
  Shannon = diversity(matrix, index = "shannon"),
  Simpson = diversity(matrix, index = "simpson"),
  Richness = rowSums(matrix > 0),  # Count of non-zero genera
  Evenness = diversity(matrix, index = "shannon") / log(rowSums(matrix > 0))  # Pielou’s Evenness
)

write_tsv(alpha_diversity, file.path(work_dir, "tables", paste0("alpha_diversities.tsv")))

# Beta Diversity Measures
bray_dist <- vegdist(matrix , method = "bray")
bray_mat <- as.matrix(bray_dist)
write.table(bray_mat, file = file.path(work_dir, "tables", paste0("bray_dist_matrix.tsv")), sep = "\t", quote = FALSE)


jaccard_dist <- vegdist(matrix , method = "jaccard", binary = TRUE)
jaccard_mat <- as.matrix(jaccard_dist)
write.table(jaccard_mat, file.path(work_dir, "tables", paste0("jaccard_dist_matrix.tsv")), sep = "\t", quote = FALSE)
