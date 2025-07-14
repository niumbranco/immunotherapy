#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

# Option parsing
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing filtered FASTQ files"),
  make_option(c("-o", "--output_dir"), type="character", default="dada2_out", help="Directory to save outputs"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads to use")
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

# Load ASV inference results (loading denoised objects)
dadaFs <- readRDS(file.path(opt$output_dir, "dadaFs.rds"))
dadaRs <- readRDS(file.path(opt$output_dir, "dadaRs.rds"))

# Retrieve sample names from filtered files (filtered FASTQ paths)
filtFs <- sort(list.files(opt$input_dir, pattern = "_F_filt.fastq.gz$", full.names = TRUE))
filtRs <- sort(list.files(opt$input_dir, pattern = "_R_filt.fastq.gz$", full.names = TRUE))

# Extract sample names 
sample.names <- gsub("_F_filt.fastq.gz", "", basename(filtFs))
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

# Merge paired reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for (i in seq_along(sample.names)) {
  sample <- sample.names[i]
  cat("Processing sample:", sample, "\n")
  
  derepF <- readRDS(file.path(opt$output_dir, "derep_objects", paste0(sample, "_F_derep.rds")))
  derepR <- readRDS(file.path(opt$output_dir, "derep_objects", paste0(sample, "_R_derep.rds")))
  
  merger <- mergePairs(dadaFs[[sample]], derepF, dadaRs[[sample]], derepR)
  mergers[[sample]] <- merger
  
  rm(derepF, derepR)
  gc()
}

# Create sequence table
seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")
saveRDS(seqtab, file.path(opt$output_dir, "seqtab.rds"))

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = opt$threads)
cat("Non-chimeric table dimensions:", dim(seqtab.nochim), "\n")
saveRDS(seqtab.nochim, file.path(opt$output_dir, "seqtab_nochim.rds"))