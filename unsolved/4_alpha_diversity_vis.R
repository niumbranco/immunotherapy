library(tidyverse)
library(vegan)
library(ape)

# working directory # where to save the results
work_dir <- "/Users/smartinezarbas/Repositories/nium/friesland_campina/data"
setwd(work_dir)

alpha_diversity <- read_tsv(file.path(work_dir, "tables", "alpha_diversities.tsv"), show_col_types = F)
metadata <- read_tsv(file.path(work_dir,"tables", "metadata.tsv"), show_col_types = FALSE)

# Add group column
alpha_diversity <- alpha_diversity %>% 
  left_join(., metadata, by = c("Sample_ID" = "sample"))

## alpha diversity trends
alpha_trends_plot <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness, Evenness), 
               names_to = "Diversity_Measure", 
               values_to = "Value") %>%
  filter(Diversity_Measure == "Shannon") %>% 
  mutate(time = as_factor(time)) %>%
  mutate(time = fct_relevel(time, c("0", "6", "12", "24"))) %>%
  mutate(condition = fct_relevel(condition, c("F", "M", "H"))) %>%
  ggplot(aes(x = time, y = Value)) +
  geom_point(width = 0.2, alpha = 0.6, color = "black") + 
  #geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "steelblue", linewidth = .7) +  # Trend line
  stat_summary(fun = mean, aes(group = 1), geom = "smooth", method = "loess", se = FALSE, color = "turquoise", linewidth = .7) +
  facet_wrap(~condition, scales = "free", ncol = 1) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Alpha Diversity by Condition (F, M, H)",
       x = "Time (h)",
       y = "Diversity Value (Shannon)")

pdf(file = file.path(work_dir, "visuals", "1_alpha_diversity_trends.pdf"),
    width = 10/2.54,
    height = 10/2.54)
print(alpha_trends_plot)
dev.off()

## richness trends
richness_trends_plot <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness, Evenness), 
               names_to = "Diversity_Measure", 
               values_to = "Value") %>%
  filter(Diversity_Measure == "Richness") %>% 
  mutate(time = as_factor(time)) %>%
  mutate(time = fct_relevel(time, c("0", "6", "12", "24"))) %>%
  mutate(condition = fct_relevel(condition, c("F", "M", "H"))) %>%
  ggplot(aes(x = time, y = Value)) +
  geom_point(width = 0.2, alpha = 0.6, color = "black") + 
  #geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "steelblue", linewidth = .7) +  # Trend line
  stat_summary(fun = mean, aes(group = 1), geom = "smooth", method = "loess", se = FALSE, color = "turquoise", linewidth = .7) +
  facet_wrap(~condition, scales = "free", ncol = 1) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Richness by Condition (F, M, H)",
       x = "Time (h)",
       y = "Richness")

pdf(file = file.path(work_dir, "visuals", "1_richness_trends.pdf"),
    width = 10/2.54,
    height = 10/2.54)
print(richness_trends_plot)
dev.off()

## alpha diversity trends
Evenness_trends_plot <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness, Evenness), 
               names_to = "Diversity_Measure", 
               values_to = "Value") %>%
  filter(Diversity_Measure == "Evenness") %>% 
  mutate(time = as_factor(time)) %>%
  mutate(time = fct_relevel(time, c("0", "6", "12", "24"))) %>%
  mutate(condition = fct_relevel(condition, c("F", "M", "H"))) %>%
  ggplot(aes(x = time, y = Value)) +
  geom_point(width = 0.2, alpha = 0.6, color = "black") + 
  #geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "steelblue", linewidth = .7) +  # Trend line
  stat_summary(fun = mean, aes(group = 1), geom = "smooth", method = "loess", se = FALSE, color = "turquoise", linewidth = .7) +
  facet_wrap(~condition, scales = "free", ncol = 1) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Evenness by Condition (F, M, H)",
       x = "Time (h)",
       y = "Evenness")

pdf(file = file.path(work_dir, "visuals", "1_evenness_trends.pdf"),
    width = 10/2.54,
    height = 10/2.54)
print(Evenness_trends_plot)
dev.off()