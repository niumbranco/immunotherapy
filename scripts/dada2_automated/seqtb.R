#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

# Option parsing
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Directory containing filtered FASTQ files"),
  make_option(c("-o", "--output_dir"), type = "character", default = "dada2_out", help = "Directory to save outputs"),
  make_option(c("-t", "--threads"), type = "integer", default = 4, help = "Number of threads to use")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load ASV inference results
dadaFs <- readRDS(file.path(opt$output_dir, "dadaFs.rds"))
dadaRs <- readRDS(file.path(opt$output_dir, "dadaRs.rds"))

# Retrieve filtered file paths
filtFs <- sort(list.files(opt$output_dir, pattern = "_16S_F_filt.fastq.gz$", full.names = TRUE))
filtRs <- sort(list.files(opt$output_dir, pattern = "_16S_R_filt.fastq.gz$", full.names = TRUE))

# Extract sample names
sample.names <- gsub("_F_filt.fastq.gz$", "", basename(filtFs))
print(sample.names)
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

# Merge paired reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for (i in seq_along(sample.names)) {
  sample <- sample.names[i]
  cat("Processing sample:", sample, "\n")
  
  derepF_path <- file.path(opt$output_dir, "derep_objects", paste0(sample, "_F_derep.rds"))
  derepR_path <- file.path(opt$output_dir, "derep_objects", paste0(sample, "_R_derep.rds"))
  
  if (!file.exists(derepF_path) || !file.exists(derepR_path)) {
    cat("⚠️ Missing derep files for sample:", sample, "\n")
    next
  }
  
  derepF <- readRDS(derepF_path)
  derepR <- readRDS(derepR_path)
  
  merger <- mergePairs(dadaFs[[sample]], derepF, dadaRs[[sample]], derepR)
  cat(" → Merged dim: ", paste(dim(merger), collapse = " × "), "\n")
  
  mergers[[sample]] <- merger
  rm(derepF, derepR)
  gc()
}

# Check mergers
non_empty_mergers <- mergers[!sapply(mergers, function(x) is.null(x) || nrow(x) == 0)]
cat("Number of non-empty mergers:", length(non_empty_mergers), "/", length(mergers), "\n")
saveRDS(mergers, file.path(opt$output_dir, "debug_mergers.rds"))

# Create sequence table
if (length(non_empty_mergers) == 0) {
  stop("All mergers are empty — cannot create sequence table.")
}

seqtab <- makeSequenceTable(non_empty_mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")
saveRDS(seqtab, file.path(opt$output_dir, "seqtab.rds"))

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = opt$threads)
cat("Non-chimeric table dimensions:", dim(seqtab.nochim), "\n")
saveRDS(seqtab.nochim, file.path(opt$output_dir, "seqtab_nochim.rds"))