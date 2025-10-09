#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

# --- Inputs and outputs from Snakemake ---
json_files <- snakemake@input
out_csv    <- snakemake@output[[1]]
out_pdf    <- snakemake@output[[2]]

# --- Collect and summarize metrics ---
results_list <- lapply(json_files, function(file) {
  sample_id <- gsub("_report.json", "", basename(file))
  j <- fromJSON(file)
  before_reads <- j$summary$before_filtering$total_reads
  after_reads  <- j$summary$after_filtering$total_reads
  q20          <- j$summary$after_filtering$q20_rate
  q30          <- j$summary$after_filtering$q30_rate
  data.frame(sample_id, reads_before = before_reads, reads_after = after_reads,
             q20_rate = q20, q30_rate = q30)
})

summary_df <- bind_rows(results_list)

# --- Quality flags ---
min_reads <- 30000
min_q20 <- 0.90
min_q30 <- 0.85

summary_df <- summary_df %>%
  mutate(status = ifelse(
    reads_after >= min_reads & q20_rate >= min_q20 & q30_rate >= min_q30,
    "ok", "low_quality"
  ),
  reads_removed = reads_before - reads_after)

# --- Write summary table ---
write.csv(summary_df, out_csv, row.names = FALSE)

# --- Plots ---
p_reads <- ggplot(summary_df, aes(x = sample_id)) +
  geom_bar(aes(y = reads_after), stat = "identity", fill = "skyblue") +
  coord_flip() + theme_bw() +
  labs(title = "Reads after filtering", y = "Reads", x = "Sample")

p_q20 <- ggplot(summary_df, aes(x = q20_rate)) +
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "darkgreen") +
  geom_vline(xintercept = min_q20, color = "red", linetype = "dashed") +
  theme_bw() + labs(title = "Q20 distribution", x = "Q20 rate")

p_q30 <- ggplot(summary_df, aes(x = q30_rate)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "blue") +
  geom_vline(xintercept = min_q30, color = "red", linetype = "dashed") +
  theme_bw() + labs(title = "Q30 distribution", x = "Q30 rate")

p_status <- ggplot(summary_df, aes(x = status, fill = status)) +
  geom_bar() + theme_bw() +
  labs(title = "Sample quality categories", x = "", y = "Count") +
  scale_fill_manual(values = c("ok" = "springgreen3", "low_quality" = "firebrick2"))

pdf(out_pdf, width = 10, height = 8)
print(ggarrange(p_reads, p_q20, p_q30, p_status, ncol = 2, nrow = 2))
dev.off()

cat("✅ fastp summary complete.\n")
cat("Samples processed:", nrow(summary_df), "\n")
cat("Low quality samples:", sum(summary_df$status == "low_quality"), "\n")
