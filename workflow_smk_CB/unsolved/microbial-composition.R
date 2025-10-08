##1.Load packages
library(phyloseq)
library(ggplot2)

##2. load the phyloseq object from dada2 analysis 
ps <- readRDS("ps.rds")

##3. Check the taxonomic ranks available 
tax_table(ps)

##4. Plot absolute abundance barplot by Phylum => shows the number of sequences per sample colored by Phylum 
plot_bar(ps, fill="Phylum") + 
  theme_minimal() +
  labs(title="Microbial Composition by Phylum", x="Samples", y='Number of sequences') +
  theme(axis.text.x=element_blank())

##5. Normalize abundances to relative proportions 
ps.prop <- transform_sample_counts(ps, function(x) x / sum(x)) #each sample is scaled to that total abundance = 1 

##6. Plot relative abundance by Phylum => represent the proportion of each phylum in ONE sample 
plot_bar(ps.prop, fill = "Phylum") +
  theme_minimal() +
  labs(title = "Relative Microbial Composition by Phylum", x = "Samples", y = "Proportion") +
  theme(axis.text.x = element_blank())


##7. Plot relative abundance by Genus 
#there are too many Genera so the plot becomes unreadable 
#so we keep only the top 10 most abundant genera across all samples
ps.genus <- tax_glom(ps.prop, taxrank="Genus")
top10_genus <- names(sort(taxa_sums(ps.genus), decreasing = TRUE))[1:10]
ps.top10 <- prune_taxa(top10_genus, ps.genus)
plot_bar(ps.top10, fill="Genus") + 
  theme_minimal() +
  labs(title="Microbial Composition by Genus", x="Samples", y='Number of sequences') +
  theme(axis.text.x=element_blank())

#plot relative abundance by Genus but without verifying that these taxa are genus one 
top10_genus <- names(sort(taxa_sums(ps.prop), decreasing = TRUE))[1:10]
ps.top10 <- prune_taxa(top10_genus, ps.prop)
plot_bar(ps.top10, fill="Genus") + 
  theme_minimal() +
  labs(title="Microbial Composition by Genus", x="Samples", y='Number of sequences') +
  theme(axis.text.x=element_blank())
