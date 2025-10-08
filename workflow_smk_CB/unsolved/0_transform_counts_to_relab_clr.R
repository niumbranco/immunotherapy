
# upload packages
library(tidyverse)
library(compositions)

# working directory # where to save the results
work_dir <- "/Users/smartinezarbas/Repositories/nium/friesland_campina/data/tables"
setwd(work_dir)

# check data: make relabs for Phylum, to visualize, make relabs for everything, and CLR, and save
ranks <- c("P", "G", "S") 

# prefiltering of data
remove_counts_of <- c("Chordata", "Homo", "Homo sapiens")

# remove also very low counts (less than 10 counts)
for (i in 1:length(ranks)){
  
  rank_i <- ranks[i]
  taxa_in <- file.path(work_dir, paste0("taxonomy_counts_", rank_i, ".tsv"))
  taxa_counts <- read_tsv(taxa_in, show_col_types = F)
  
  taxa_filt <- taxa_counts %>% 
    filter(!name %in% remove_counts_of)
  
  taxa_filt %>% 
    map_if(is_double, function(x) ifelse(is.na(x), 0, x)) %>% as_tibble() %>% 
    map_if(is_double, function(x) ifelse(x <= 10, 0, x)) %>% as_tibble() %>% # remove counts < 10
    gather(sample, counts, -name) %>% 
    group_by(sample) %>% 
    mutate(relab = counts * 100 / sum(counts)) %>% 
    select(-counts) %>% 
    spread(sample, relab) %>% 
    write_tsv(., file.path(work_dir, paste0("Relabs_", rank_i, ".tsv")))

  # transform to CLR and save
  
  taxa_filt %>% 
    map_if(is_double, function(x) ifelse(is.na(x), 0, x)) %>% as_tibble() %>% 
    gather(sample, counts, -name) %>% 
    group_by(sample) %>% 
    mutate(relab = counts * 100 / sum(counts)) %>% 
    select(-counts) %>% 
    spread(sample, relab) %>% 
    column_to_rownames("name") %>% 
    #as.matrix() %>% 
    compositions::clr(.) %>% 
    as_tibble() %>% 
    mutate(name = taxa_filt$name, .before = 1) %>% 
    write_tsv(., file.path(work_dir, paste0("CLR_", rank_i, ".tsv")))
  
}
