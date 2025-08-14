#!/usr/bin/env Rscript

#load libraries 
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)

##------------to execute the script alone (without snakemake)--------------
seqtab_file <- "~/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/PRJEB61942/reads_taxonomy/seqtab.rds"
taxa_file <- "~/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/PRJEB61942/reads_taxonomy/taxa_fixed.rds"
metadata_file <- "~/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/metadata/SraRunTable.csv"
##-------------------------------------------------------------------------

path <- "~/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/fastp_output"

#read arguments
args <- commandArgs(trailingOnly = TRUE)
#seqtab_file <- args[1]    #seqtab.rds
#taxa_file <- args[2]      #taxa.rds
#metadata_file <- args[3]  #sraRunTable.csv
alpha_out <- args[4]      #output alpha diversity file
beta_out <- args[5]       #output beta diversity file

alpha_out <- "alpha_div.tsv"
beta_out <- "bray_matrix.tsv"

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

# Attribution des noms de lignes à la table de séquences
#rownames(seqtab.nochim) <- sample.names

# Attribution des noms de lignes à la table de séquences
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
  warning("Des ASV dupliqués ont été trouvés dans taxa. On garde la première occurrence.")
  taxa_raw <- taxa_raw[!duplicated(rownames(taxa_raw)), ]
}

# Conserver uniquement les ASV communs
common_asv <- intersect(colnames(seqtab.nochim), rownames(taxa_raw))
seqtab.nochim <- seqtab.nochim[, common_asv]
taxa_raw <- taxa_raw[common_asv, ]

# Créer un objet tax_table compatible avec phyloseq
taxa <- tax_table(as.matrix(taxa_raw))

cat("Assay.Type uniques :", unique(metadata$Assay.Type), "\n")

#prepare metadata
##metadata.16S <- metadata[metadata$Assay.Type == "AMPLICON", ]
metadata.16S <- metadata[grepl("AMPLICON", metadata$Assay.Type, ignore.case = TRUE), ]
metadata.16S$SampleID <- metadata.16S$Run
metadata.16S <- metadata.16S[metadata.16S$SampleID %in% sample.names, ]

samdf <- data.frame(
  row.names = metadata.16S$SampleID,
  SampleID = metadata.16S$SampleID,
  #Group = sub("-.*", "", metadata.16S$Submitter_Id)
  Group = "AMPLICON"
)

common_samples <- intersect(rownames(seqtab.nochim), rownames(samdf))

seqtab.nochim <- seqtab.nochim[common_samples, ]
samdf <- samdf[common_samples, ]

#build phyloseq object 
#rownames(seqtab.nochim) <- sample.names
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

### PLOT ALPHA DIVERSITY 
library(ggplot2)
alpha_div <- read.table("alpha_div.tsv", header = TRUE, sep = "\t")

# SHANNON 
ggplot(alpha_div, aes(x = Group, y = Shannon)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.1, color = "darkblue") +
  labs(title = "Diversité alpha (Shannon)", x = "Groupe", y = "Indice de Shannon") +
  theme_minimal()

# SIMPSON
ggplot(alpha_div, aes(x = Group, y = Simpson)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.1, color = "darkblue") +
  labs(title = "Diversité alpha (Simpson)", x = "Groupe", y = "Indice de Simpson") +
  theme_minimal()

### PLOT BETA DIVERISTY

# Lire correctement la matrice
bray <- read.table("bray_matrix.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
bray <- as.matrix(bray)

# Forcer les noms des lignes à être des chaînes
rownames(bray) <- as.character(rownames(bray))

# Corriger la matrice pour que noms lignes = colonnes
common_ids <- intersect(rownames(bray), colnames(bray))
bray <- bray[common_ids, common_ids]

# Nettoyer les éventuels NA
bray[is.na(bray)] <- 0

# Vérifier avant PCoA
print(dim(bray))  # doit afficher (87, 87)

# Créer la distance et faire la PCoA
bray_dist <- vegdist(bray, method = "bray")
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Construire un dataframe pour ggplot
df_pcoa <- data.frame(
  SampleID = rownames(pcoa$points),
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2]
)

# Ajouter les groupes
alpha_div <- read.table("alpha_div.tsv", header = TRUE, sep = "\t")
df_pcoa$SampleID <- as.character(df_pcoa$SampleID)
alpha_div$SampleID <- as.character(alpha_div$SampleID)

df_pcoa <- merge(df_pcoa, alpha_div[, c("SampleID", "Group")], by = "SampleID")

# Tracer le PCoA
print(
  ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3) +
    labs(title = "PCoA (Bray-Curtis)", x = "PC1", y = "PC2") +
    theme_minimal()
)

# Barplot: microbial composition by Phylum
ps.phylum <- tax_glom(ps.prop, taxrank = "Phylum")
phylum_df <- psmelt(ps.phylum) #transformer en dataframe pour ggplot
phylum_df$Phylum <- as.character(phylum_df$Phylum)

# Barplot: microbial composition by Genus
ps.genus <- tax_glom(ps.prop, taxrank = "Genus")
genus_df <- psmelt(ps.genus)
rare_genus <- genus_df %>%  # Regrouper les genres très peu abondants
  group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  filter(total < 0.01) %>%   
  pull(Genus)
genus_df$Genus[genus_df$Genus %in% rare_genus] <- "Other"

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
p_pcoa <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCoA (Bray-Curtis)", x = "PC1", y = "PC2") +
  theme_classic()

# Microbial composition by Phylum 
p_phylum <- ggplot(phylum_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative abundance by Phylum", y = "Abundance", x = NULL) +
  theme_classic() +
  theme(axis.text.X = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette="Set3")

# Microbial composition by Genus
nb_colors <- length(unique(genus_df$Genus))
custom_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(nb_colors)
genus_df$Genus <- droplevels(factor(genus_df$Genus))
p_genus <- ggplot(genus_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative abundance by Genus", y = "Abundance", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 6),       # pour éviter que la légende déborde
    legend.key.size = unit(0.4, "cm"),
  ) +
  scale_fill_manual(values = custom_palette)

# COMBINER
combined_plot <- ggarrange(
  ggarrange(p_shannon, p_simpson, ncol = 2, labels = c("A", "B")),
  p_pcoa,
  p_phylum,
  p_genus,
  nrow = 4,
  labels = c("", "C", "D", "E")
)

print(combined_plot)

# Sauvegarde avec fond blanc
ggsave(
  filename = "/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/PRJEB61942/results/diversity_combined_plot_all_with_Genus.pdf", 
  plot= combined_plot, 
  width = 18, #largeur
  height = 23, #hauteur
  bg = "white",
  units = "in",
  dpi = 300
  )
