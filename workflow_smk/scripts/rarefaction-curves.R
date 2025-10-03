#################################

##RAREFACTION CURVES

library(phyloseq)
library(parallel)
library(doParallel)
library(vegan)

ps <- readRDS("ps.rds") #load phyloseq object 

##1. Define the maximum sequencing depth across all samples 
max.depth <- max(sample_sums(ps))
max.depth   #gives the maximum number of reads in any sample 

##2. Create a sequence of rarefaction depths (here 250 steps between 2 and max.depth)
sampling.depth <- round(seq(from=2, to=max.depth, length.out=250))

##3. creation of the OTU table 
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) {
  otu <- t(otu)
}

##4. set up parallel processing  
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

##5. compute rarefaction data for each sample
rarecurve_data <- foreach(i=1:nrow(otu), .combine=rbind) %dopar% {
  row <- otu[i, ]
  rarefactions <- sapply(sampling.depth, function(d) {
    if (sum(row) < d) {
      return(NA)
    }
    else{
      return(vegan::rarefy(row, sample=d))
    }
  })
  data.frame(SampleID = rownames(otu)[i],
             Depth=sampling.depth, 
             Richness=rarefactions)
}

stopCluster(cl)

##6. Plot rarefaction curves colored by group
library(ggplot2)
library(reshape2)

cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

##TO COLOR BY GROUP: 

#prepare metadata for merging 
sample_info <- data.frame(sample_data(ps))
sample_info$SampleID <- rownames(sample_info)
#merge rarefaction data with metadata
rarecurve_data_merged <- merge(rarecurve_data_clean, sample_info, by ="SampleID") 
#plot 
ggplot(rarecurve_data_merged, aes(x=Depth, y=Richness, group=SampleID, color=Group)) +
  geom_line(alpha=0.5) +
  labs(title="Rarefaction curves by group", x="Sequencing depth", y="species richness") + 
  theme_minimal()

#additional filtering for comparison
otu_raw <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_raw <- t(otu_raw) }

taxa_singletons <- rowSums(otu_raw > 0) == 1
otu_no_singletons <- otu_raw[!taxa_singletons, ]

taxa_single_reads <- rowSums(otu_raw) == 1 
otu_no_single_reads <- otu_raw[!taxa_single_reads, ]

#function to compute rarefaction data
compute_rarefaction_data <- function(otu_table, sampling.depth) {
  foreach(i=1:nrow(otu_table), .combine=rbind) %dopar% {
    row <- otu_table[i, ]
    rarefactions <- sapply(sampling.depth, function(d) {
      if (sum(row) < d) return (NA) 
      else return(vegan::rarefy(row, sample=d))
    })
    data.frame(SampleID = rownames(otu_table)[i], 
               Depth=sampling.depth, 
               richness=rarefactions)
  }
}

#apply rarefaction to each version of the OTU table 
rare_raw <- compute_rarefaction_data(otu_raw, sampling.depth)
rare_nosingleton <- compute_rarefaction_data(otu_no_singletons, sampling.depth)
rare_nosingleread <- compute_rarefaction_data(otu_no_single_reads, sampling.depth)

stopCluster(cl)

##READ NORMALIZATION

ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)) #normalize read counts by total sum per sample (relative abundance)
ps.prop

##READ TRANSFORMATION (CLR)

#install.packages("zCompositions")
#install.packages("compositions")
library(zCompositions)
library(compositions)

#filter taxonomy table to remove unclassified phyla
tax_table <- data.frame(tax_table(ps))
cleaned_tax_table <- tax_table[!grepl("unclassified", tax_table$Phylum, ignore.case=TRUE), ]
#reconstruct phyloseq object with filtered taxonomy 
ps.filtered <- merge_phyloseq(otu_table(ps), sample_data(ps), tax_table(as.matrix(cleaned_tax_table)))
#agglomerate OTUs at the genus level
ps.genus <- tax_glom(ps.filtered, taxrank = "Genus")
ps.genus
#extract the OTU table from the phyloseq object
otu.genus <- as.data.frame(otu_table(ps.genus))
if (taxa_are_rows(ps.genus)) {
  otu.genus <- t(otu.genus)
}
#define the CLR transformation function 
transform.clr <- function(otu.df){
  print(dim(otu.df))
  d.1 <- data.frame(otu.df[which(rowSums(otu.df) > 10), ], check.names=F)
  print(dim(d.1))
  d.czm <- t(cmultRepl(t(d.1), label=0, method="CZM"))
  print(dim(d.czm))
  d.clr <- apply(d.czm, 2, function(x){log(x) - mean(log(x))})
  print(dim(d.clr))
  return(d.clr)
}
#apply CLR transformation  
gen.clr <- transform.clr(otu.genus)
gen.clr

## PCA ON CLR TRANSFORMED GENUS LEVEL TABLE 

#keep only samples present in phyloseq object 
valid_samples <- sample_names(ps)
gen.clr.filtered <- gen.clr[rownames(gen.clr) %in% valid_samples, ]

#remove constant columns 
gen.clr.clean <- gen.clr.filtered[, apply(gen.clr.filtered, 2, function(x) sd(x) > 0)]

#perform PCA 
pca <- prcomp(gen.clr.clean, center = TRUE, scale. = TRUE)

#create dataframe for plotting with metadata
sample_info <- data.frame(sample_data(ps))
sample_info$SampleID <- rownames(sample_info)
pca_df <- data.frame(SampleID = rownames(pca$x), pca$x)
pca_df <- merge(pca_df, sample_info[, c("SampleID", "Group")], by = "SampleID")

#Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA on CLR-transformed genus level abundance",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()

## NMDS WITH CLR TRANSFORMATION AND COLOR BY GROUP

dist.clr <- dist(gen.clr.clean, method="euclidean")

nmds.clr <- metaMDS(dist.clr, k=2, trymax=100)

nmds_scores <- as.data.frame(scores(nmds.clr))
nmds_scores$SampleID <- rownames(nmds_scores)

sample_info <- data.frame(sample_data(ps))
sample_info$SampleID <- rownames(sample_info)

nmds_df <- merge(nmds_scores, sample_info[, c("SampleID", "Group")], by="SampleID")

ggplot(nmds_df, aes(x=NMDS1, y=NMDS2, color=Group)) + 
  geom_point(size=3) +
  labs(title="NMDS on CLR transformed genus abundances (Euclidean distance)") +
  theme_minimal()
