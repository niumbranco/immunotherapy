library(dada2)
packageVersion("dada2")

##1. GETTING READY

#définir le chemin vers les fastq nettoyés avec fastp
path <- "/Users/Emma/Documents/EIDD/2A/STAGE 2A/DOCUMENTS DE STAGE /CODES/fastp_output"
list.files(path)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path, pattern="_1.clean.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.clean.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- paste0(sample.names, "_16S")

head(fnFs)
head(fnRs)
head(sample.names)

##2. INSPECT READ QUALITY PROFILES
#quality of the forward reads
plotQualityProfile(fnFs[1:3])
#quality of the reverse reads
plotQualityProfile(fnRs[1:3])

##3. LEARN THE ERROR RATES      #take time 
if (!file.exists("errF.rds") || !file.exists("errR.rds")) {   #to not recalculate if the file already exist 
  errF <- learnErrors(fnFs, multithread = TRUE)
  errR <- learnErrors(fnRs, multithread = TRUE)
  saveRDS(errF, file="errF.rds")
  saveRDS(errR, file="errR.rds")
} else {
  errF <- readRDS("errF.rds")
  errR <- readRDS("errR.rds")
}
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##4. SAMPLE INFERENCE     #take time 
if (!file.exists("dadaFs.rds") || !file.exists("dadaRs.rds")) {
  dadaFs <- dada(fnFs, err=errF, multithread = TRUE)
  dadaRs <- dada(fnRs, err=errR, multithread = TRUE)
  saveRDS(dadaFs, "dadaFs.rds")
  saveRDS(dadaRs, "dadaRs.rds")
} else {
  dadaFs <- readRDS("dadaFs.rds")
  dadaRs <- readRDS("dadaRs.rds")
}

##5. MERGE PAIRED READS
mergers <- mergePairs(dadaFs, fnFs, dadaRs, fnRs, verbose=TRUE)
head(mergers[[1]])

##6. CONSTRUCT SEQUENCE TABLE 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

##7. REMOVE CHIMERAS.     #take time 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

##8. TRACK READS THROUGH THE PIPELINE 
out <- data.frame(row.names=sample.names, input=rep(NA, length(sample.names)), filtered=rep(NA, length(sample.names)))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##9. ASSIGN TAXONOMY
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/EIDD/2A/STAGE 2A/DOCUMENTS DE STAGE /silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

##10. EVALUATE ACCURACY
#we don't need to do this step because our dataset doesn't contain a mock community

##11. BONUS: HANDOFF TO PHYLOSEQ 
library(phyloseq)
library(Biostrings)
library(ggplot2)

theme_set(theme_bw())

sample.names <- gsub("_1\\.clean\\.fastq\\.gz", "_16S", rownames(seqtab.nochim))

metadata_path <- "/Users/Emma/Documents/EIDD/2A/STAGE 2A/DOCUMENTS DE STAGE /SraRunTable.csv"
metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

metadata.16S <- metadata[metadata$Assay.Type=="AMPLICON",]
metadata.16S$SampleID <- paste0(metadata.16S$Run, "_16S")
metadata.16S <- metadata.16S[metadata.16S$SampleID %in% sample.names,]

samdf <- data.frame(
  row.names = metadata.16S$SampleID, 
  SampleID = metadata.16S$SampleID, 
  Group = sub("-.*", "", metadata.16S$Submitter_Id)
)

rownames(seqtab.nochim) <- sample.names

identical(sort(rownames(seqtab.nochim)), sort(rownames(samdf)))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa))

dna <- Biostrings::DNAStringSet((taxa_names(ps)))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, file="ps.rds")

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

#visualize alpha diversity
plot_bar(ps.prop, fill="Phylum")
plot_bar(ps.prop, fill="Phylum", x="SampleID") +
  facet_wrap(~ Group, scales = "free_x") + 
  theme(axis.text.x = element_blank())

plot_heatmap(transform_sample_counts(ps, function(x) log10(x+1)))

plot_richness(ps.prop, measures = c("Shannon", "Simpson")) 
plot_richness(ps.prop, x="Group", measures=c("Shannon", "Simpson")) +
  geom_boxplot(aes(fill=Group), alpha=0.4) + 
  labs(title="Alpha diversity (relative abundance)", y="Diversity Index", x="Group")

#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Group", title="Bray NMDS") +
  geom_point(size=3) +
  labs(title = "NMDS (Bray-Curtis)", color="Group")

#bar plot
top20 <- names(sort(taxa_sums(ps.prop), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.prop)
plot_bar(ps.top20, x="SampleID", fill="Family") +
  labs(title="Top 20 taxa (relative abundance)", y="Proportion", x="Sample") 
