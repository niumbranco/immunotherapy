#!/usr/bin/env Rscript

# Load packages
suppressPackageStartupMessages({
  library(optparse)
  library(dada2)
})

# Define options
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing trimmed FASTQ files"),
  make_option(c("-o", "--output_dir"), type="character", default="dada2_out", help="Output directory [default %default]"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of threads [default %default]"),
)

# Parse args
opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory if not exists
dir.create(file.path(opt$output_dir, "derep_objects"), showWarnings = FALSE)

# Log version
cat("DADA2 version:", as.character(packageVersion("dada2")), "\n")

# Load filtered FASTQ file paths (same logic as in err.R)
filtFs <- sort(list.files(opt$input_dir, pattern = "_F_filt.fastq.gz$", full.names = TRUE))
filtRs <- sort(list.files(opt$input_dir, pattern = "_R_filt.fastq.gz$", full.names = TRUE))
sample.names <- gsub("_F_filt.fastq.gz", "", basename(filtFs))

# Define paths to error models previously learned 
errF_path <- file.path(opt$output_dir, "errF.rds")
errR_path <- file.path(opt$output_dir, "errR.rds")
errF <- readRDS(errF_path)
errR <- readRDS(errR_path)

# initialize output lists 
dadaFs <- vector("list", length(sample.names))
dadaRs <- vector("list", length(sample.names))
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

# denoising loop 
for (i in seq_along(sample.names)) {
  sample <- sample.names[i]
  cat("Processing sample:", sample, "\n")
  
  derepF <- dada2::derepFastq(filtFs[i])
  derepR <- dada2::derepFastq(filtRs[i])
  
  saveRDS(derepF, file = file.path(opt$output_dir, "derep_objects", paste0(sample, "_F_derep.rds")))
  saveRDS(derepR, file = file.path(opt$output_dir, "derep_objects", paste0(sample, "_R_derep.rds")))
  
  dadaFs[[i]] <- dada(derepF, err = errF, multithread = opt$threads)
  dadaRs[[i]] <- dada(derepR, err = errR, multithread = opt$threads)
  
  rm(derepF, derepR)
  gc()
}

# save denoising results 
saveRDS(dadaFs, file.path(opt$output_dir, "dadaFs.rds"))
saveRDS(dadaRs, file.path(opt$output_dir, "dadaRs.rds"))