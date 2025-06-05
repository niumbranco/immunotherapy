#install.packages("jsonlite")
#install.packages("dplyr")

library(jsonlite)   #used to read and write JSON files 
library(dplyr)      #used for data manipulation like filtering, selecting, summarizing...) 

setwd("/Users/Emma/Documents/EIDD/2A/STAGE 2A/DOCUMENTS DE STAGE /CODES")

json_files <- list.files ("fastp_reports", pattern="_report.json$", full.names=TRUE) #list all JSON files in the fastp_reports directory

results_list <- list() #initialize empty list to collect results

#loop through each json files
for (file in json_files) {
  sample_id <- gsub("_report.json", "", basename(file)) #extract sample ID from the file name (ex:ERR13827976_report.json -> ERR13827976)
  json_data <- fromJSON(file) #read JSON content
  
  #extract important parameters 
  before_reads <- json_data$summary$before_filtering$total_reads
  after_reads <- json_data$summary$after_filtering$total_reads
  q20 <- json_data$summary$after_filtering$q20_rate
  q30 <- json_data$summary$after_filtering$q30_rate
  
  #create a row in a table with all these values
  row <- data.frame(sample_id=sample_id, reads_before=before_reads, reads_after=after_reads, q20_rate=q20, q30_rate=q30)
  
  #add row to the list 
  results_list[[sample_id]] <- row
}

summary_df <- bind_rows(results_list) #combine all rows into one data frame

#quality check 
min_reads <- 30000
min_q20 <- 0.90
min_q30 <- 0.85

summary_df <- summary_df %>%
  mutate(status=ifelse(
    reads_before >= min_reads & reads_after >= min_reads & q20_rate >= min_q20 & q30_rate >= min_q30, "Ok","low quality"
    )
  )

### CSV TABLE: 

write.csv(summary_df, "fastp_summary_table.csv", row.names=FALSE) #save to csv

### PRINTING INFORMATION TO THE CONSOLE 

print(summary_df) 
cat("total samples: ", nrow(summary_df), "\n")
cat("samples which not respect criteria: ", sum(summary_df$status == "low quality"), "\n")

### VISUALIZATION

library(ggplot2)

#histogram of reads_before 
ggplot(data=summary_df, aes(x=reads_before)) + 
  geom_histogram(binwidth=4000, fill="pink", color="hotpink") +
  geom_vline(xintercept=30000, color="blue", linetype="dashed") +
  annotate("text", x=30000, y=6, label= " Minimum threshold", color="blue", hjust=0) +
  labs(title="Reads before filtering", x="reads_before", y="Number of samples") +
  theme_bw()

#histogram of reads_after
ggplot(data=summary_df, aes(x=reads_after)) + 
  geom_histogram(binwidth=4000, fill="pink", color="hotpink") +
  geom_vline(xintercept=30000, color="blue", linetype="dashed") +
  annotate("text", x=30000, y=6, label= " Minimum threshold", color="blue", hjust=0) +
  labs(title="Reads after filtering", x="reads_after", y="Number of samples") +
  theme_bw()

#histogram of q20 rate 
ggplot(data=summary_df, aes(x=q20_rate)) + 
  geom_histogram(binwidth=0.01, fill="skyblue", color="steelblue") +
  geom_vline(xintercept =0.90, color="tomato", linetype="dashed") +
  labs(title="Distribution of Q20 rate", x="Q20 rate", y="Number of samples") +
  theme_bw()

#histogram of q30 rate
ggplot(data=summary_df, aes(x=q30_rate)) + 
  geom_histogram(binwidth=0.01, fill="skyblue", color="steelblue") +
  geom_vline(xintercept =0.85, color="tomato", linetype="dashed") +
  labs(title="Distribution of Q30 rate", x="Q30 rate", y="Number of samples") +
  theme_bw()

#histogram of q20 rate with density
ggplot(data=summary_df, aes(x=q20_rate)) +
  geom_density(fill="lavender", color="mediumpurple3") +
  geom_vline(xintercept=0.9, color="tomato1", linetype="dashed") +
  labs(title="Distribution of Q20 rate", x="Q20 rate", y="Number of samples") +
  theme_bw()

#histogram of q30 rate with density
ggplot(data=summary_df, aes(x=q30_rate)) +
  geom_density(fill="lavender", color="mediumpurple3") +
  geom_vline(xintercept=0.85, color="tomato1", linetype="dashed") +
  labs(title="Distribution of Q30 rate", x="Q30 rate", y="Number of samples") +
  theme_bw()

#boxplot of q20 rate
ggplot(data=summary_df, aes(x="Q20", y=q20_rate)) + 
  geom_boxplot(fill="darkolivegreen2") +
  labs(title="Boxplot of Q20 rate", x="", y="Q20 rate") +
  theme_bw()

#boxplot of q30 rate
ggplot(data=summary_df, aes(x="Q30", y=q30_rate)) + 
  geom_boxplot(fill="darkolivegreen2") +
  labs(title="Boxplot of Q30 rate", x="", y="Q30 rate") +
  theme_bw() 

#barplot with the number of samples that are "ok" vs "low quality" 
ggplot(data=summary_df, aes(x=status, fill=status)) +
  geom_bar() +
  scale_fill_manual(values=c("Ok"="springgreen3", "low quality"="firebrick2")) +
  geom_text(stat="count", aes(label=after_stat(count)), vjust=2, color="black") +
  labs(title="Sample quality barplot", x="Status", y="count") +
  theme_bw()

#boxplot of Q20 rate by their status ("ok" vs "low quality")
ggplot(summary_df, aes(x=status, y=q20_rate)) +
  geom_boxplot(fill="darkolivegreen3") +
  labs(title="Q20 rate by status", x="Sample quality", y="Q20 rate") +
  theme_bw()

#boxplot of Q30 rate by their status ("ok" vs "low quality")
ggplot(summary_df, aes(x=status, y=q30_rate)) +
  geom_boxplot(fill="darkolivegreen3") +
  labs(title="Q30 rate by status", x="Sample quality", y="Q30 rate") +
  theme_bw()

#scatter plot of reads_before vs reads_after after filtering using fastp
ggplot(summary_df, aes(x=reads_before, y=reads_after, color=status)) + 
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Ok"="springgreen3", "low quality"="firebrick2")) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
  labs(title="reads_before VS reads_after after filtering with fastp", x="reads before", y="reads after", color="sample status") +
  coord_cartesian(xlim = c(0,380000), ylim=c(0,380000)) +
  theme_bw()

#histogram of the reads that have been removed by fastp 
summary_df <- summary_df %>%
  mutate(reads_removed = reads_before - reads_after)

ggplot(summary_df, aes(x=reads_removed)) +
  geom_histogram(binwidth=1000, fill="pink", color="maroon") +
  labs(title="Reads removed by fastp", x="Reads removed", y="Samples") +
  theme_bw()

#print the lowest Q20 and Q30 scores and their sample IDs 
lowest_q20 <- summary_df %>% arrange(q20_rate) %>% head(1)
lowest_q30 <- summary_df %>% arrange(q30_rate) %>% head(1)

cat("Lowest Q20 scores: \n")
print(lowest_q20[, c("sample_id", "q20_rate")])

cat("Lowest Q30 scores: \n")
print(lowest_q30[, c("sample_id", "q30_rate")])

summary_df <- summary_df %>% 
  mutate(final_status = ifelse(reads_after >= 30000, "good quality", "low quality"))

table(summary_df$final_status)

ggplot(summary_df, aes(x=final_status, fill=final_status)) +
  geom_bar() +
  scale_fill_manual(values=c("good quality"="springgreen3", "low quality"="firebrick2")) +
  geom_text(stat="count", aes(label=after_stat(count)), vjust=2, color="black") +
  labs(title="Sample classification based on reads only", x="Read quality", y="Number of samples") +
  theme_bw()

#list of samples that are low quality 
({
  cat("List of samples that are low quality based on the number of reads only:\n")
  summary_df[summary_df$reads_after < 30000, c("sample_id", "reads_after")]
})
