suppressPackageStartupMessages({
  library(jsonlite)
  library(dada2)
  library(tidyverse)
})

# Input and output from Snakemake
fastp_reports <- snakemake@input[["fastp_reports"]]
dadaFs <- snakemake@input[["dadaFs"]]
dadaRs <- snakemake@input[["dadaRs"]]
mergers <- snakemake@input[["mergers"]]
seqtab_nochim <- snakemake@input[["seqtab_nochim"]]
output_tsv <- snakemake@output[["track_table"]]
output_warn <- snakemake@output[["warnings"]]

# ---- Get raw and filtered reads from fastp JSONs ----
raw_counts <- map_dfr(fastp_reports, function(fp_json) {
  j <- fromJSON(fp_json)
  tibble(
    sample = sub("_report\\.json$", "", basename(fp_json)),
    raw = j$summary$before_filtering$total_reads,
    filtered = j$summary$after_filtering$total_reads
  )
})

# ---- Extract counts from DADA2 RDS files ----
# Helper to count reads from dada2 objects
# --- Flexible helper functions for counting reads ---

# Safe version of getN that works for many object types
safe_getN <- function(obj) {
  if (inherits(obj, "dada")) {
    return(sum(getUniques(obj)))
  } else if (inherits(obj, "derep")) {
    return(sum(getUniques(obj)))
  } else if (inherits(obj, "data.frame") && all(c("sequence", "abundance") %in% colnames(obj))) {
    return(sum(obj$abundance))
  } else if (is.matrix(obj)) {
    return(sum(obj))
  } else if (is.list(obj)) {
    # Recursively handle lists of DADA2 objects
    return(sum(sapply(obj, safe_getN)))
  } else {
    return(0)
  }
}

# Wrapper that reads the .rds file and counts reads robustly
get_nreads <- function(rds_path) {
  obj <- readRDS(rds_path)
  safe_getN(obj)
}

dadaF_counts <- tibble(sample = sub("_dadaFs\\.rds$", "", basename(dadaFs)),
                       denoisedF = sapply(dadaFs, get_nreads))
dadaR_counts <- tibble(sample = sub("_dadaRs\\.rds$", "", basename(dadaRs)),
                       denoisedR = sapply(dadaRs, get_nreads))
merger_counts <- tibble(sample = sub("_mergers\\.rds$", "", basename(mergers)),
                        merged = sapply(mergers, get_nreads))

# ---- Non-chimeric counts ----
seqtab_nochim_obj <- readRDS(seqtab_nochim)

# Check if the seqtab has rownames
if (!is.null(rownames(seqtab_nochim_obj)) && length(rownames(seqtab_nochim_obj)) > 0) {
  sample_names <- rownames(seqtab_nochim_obj)
} else {
  # fallback if rownames missing — assign sample indices
  sample_names <- paste0("Sample_", seq_len(nrow(seqtab_nochim_obj)))
  warning("seqtab_nochim.rds has no rownames — using generic sample names.")
}

nonchim_counts <- tibble(
  sample = sample_names,
  nonchim = rowSums(seqtab_nochim_obj)
)


# ---- Combine all ----
track <- raw_counts %>%
  left_join(dadaF_counts, by="sample") %>%
  left_join(dadaR_counts, by="sample") %>%
  left_join(merger_counts, by="sample") %>%
  left_join(nonchim_counts, by="sample") %>%
  mutate(retained_pct = round((nonchim / raw) * 100, 2))

# ---- Save summary ----
write_tsv(track, output_tsv)

# ---- Generate warnings ----
warnings <- track %>%
  filter(retained_pct < 50) %>%
  mutate(message = paste0("⚠️  Sample ", sample, " retained only ", retained_pct, "% of reads."))

if (nrow(warnings) > 0) {
  write_lines(warnings$message, output_warn)
} else {
  write_lines("All samples retained >50% of reads.", output_warn)
}
