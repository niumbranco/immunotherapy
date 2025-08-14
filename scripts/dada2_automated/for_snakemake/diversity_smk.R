#!/usr/bin/env Rscript

#load libraries 
suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggpubr)
})

#read arguments
args <- commandArgs(trailingOnly = TRUE)
seqtab_file <- args[1]    #seqtab.rds
taxa_file <- args[2]      #taxa.rds
metadata_file <- args[3]  #sraRunTable.csv
alpha_out <- args[4]      #output alpha diversity file
beta_out <- args[5]       #output beta diversity file

path <- "~/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/fastp_output"

#load data
seqtab.nochim <- readRDS(seqtab_file)
taxa_raw <- readRDS(taxa_file)
metadata <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

raw_files <- list.files(path, pattern = "\\.fastq\\.gz$", full.names = FALSE)
sample.names <- unique(sub("_[12]\\.clean\\.fastq\\.gz", "", raw_files))

# Vérification de la correspondance
if (length(sample.names) < nrow(seqtab.nochim)) {
  warning(paste0("Seulement ", length(sample.names), " noms extraits pour ", nrow(seqtab.nochim), " lignes dans seqtab.nochim. Ajustement forcé."))
  sample.names <- sample.names[1:nrow(seqtab.nochim)]
} else if (length(sample.names) > nrow(seqtab.nochim)) {
  warning("Trop de noms extraits. On tronque.")
  sample.names <- sample.names[1:nrow(seqtab.nochim)]
}
rownames(seqtab.nochim) <- sample.names

# Corriger la structure de taxa si nécessaire
if (all(rownames(taxa_raw) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) {
  cat("Transposition forcée car les rangs sont en rownames au lieu des colonnes.\n")
  taxa_raw <- t(taxa_raw)
  colnames(taxa_raw) <- rownames(readRDS(taxa_file))  # noms des rangs
  rownames(taxa_raw) <- colnames(readRDS(taxa_file))  # noms des ASV
}

# Supprimer les ASV dupliqués dans taxa
if (any(duplicated(rownames(taxa_raw)))) {
  taxa_raw <- taxa_raw[!duplicated(rownames(taxa_raw)), ]
}

# Conserver uniquement les ASV communs
common_asv <- intersect(colnames(seqtab.nochim), rownames(taxa_raw))
seqtab.nochim <- seqtab.nochim[, common_asv]
taxa_raw <- taxa_raw[common_asv, ]

# Créer un objet tax_table compatible avec phyloseq
taxa <- tax_table(as.matrix(taxa_raw))

#prepare metadata
metadata.16S <- metadata[grepl("AMPLICON", metadata$Assay.Type, ignore.case = TRUE), ]
metadata.16S$SampleID <- metadata.16S$Run
metadata.16S <- metadata.16S[metadata.16S$SampleID %in% sample.names, ]

samdf <- data.frame(
  row.names = metadata.16S$SampleID,
  SampleID = metadata.16S$SampleID,
  Group = "AMPLICON"
)

common_samples <- intersect(rownames(seqtab.nochim), rownames(samdf))
seqtab.nochim <- seqtab.nochim[common_samples, ]
samdf <- samdf[common_samples, ]

#build phyloseq object 
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
  sample_data(samdf), 
  taxa
)

#normalize to relative abundance 
ps.prop <- transform_sample_counts(ps, function(x) x / sum(x))

#alpha diversity 
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div$Group <- sample_data(ps.prop)$Group

#beta diversity (Bray-Curtis distance matrix)
bray <- phyloseq::distance(ps.prop, method="bray")
bray_matrix <- as.matrix(bray)
sample_ids <- rownames(alpha_div)
rownames(bray_matrix) <- sample_ids
colnames(bray_matrix) <- sample_ids

#save results 
write.table(alpha_div, file = alpha_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bray_matrix, file = beta_out, sep = "\t", quote = FALSE, col.names = NA)


### PDF REPORT WITH THE VISUALISATIONS

# SHANNON
p_shannon <- ggplot(alpha_div, aes(x = Group, y = Shannon)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.1, color = "darkblue") +
  labs(title = "Shannon index", x = "Group", y = "Shannon") +
  theme_classic()

# SIMPSON
p_simpson <- ggplot(alpha_div, aes(x = Group, y = Simpson)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.1, color = "darkblue") +
  labs(title = "Simpson index", x = "Group", y = "Simpson") +
  theme_classic()

# PCOA
bray_dist <- vegdist(bray_matrix, method = "bray")
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

df_pcoa <- data.frame(
  SampleID = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

df_pcoa <- merge(df_pcoa, alpha_div[, c("SampleID", "Group")], by = "SampleID")

p_pcoa <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCoA (Bray-Curtis)", x = "PC1", y = "PC2") +
  theme_classic()

# COMBINER
combined_plot <- ggarrange(
  ggarrange(p_shannon, p_simpson, ncol = 2, labels = c("A", "B")),
  p_pcoa,
  nrow = 2,
  labels = c("", "C")
)

# Sauvegarde
ggsave(filename = sub("bray_matrix.tsv", "diversity_combined_plot.png", beta_out),
       plot = combined_plot,
       width = 10, height = 8, bg = "white")