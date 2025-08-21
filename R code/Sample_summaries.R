# Calculate means, median 
rm(list=ls())

setwd("~/MRes Biosystematics/Project 2/data/results")
getwd()

# load packages
library(tidyverse)
library(dplyr)
library(ggplot2)


sourmash_result <- read.csv("SRA_search_results.csv")

sourmash_result_pct <- sourmash_result %>% 
  mutate(p_match = f_match * 100)

total_summary_stat <- sourmash_result_pct %>%
  group_by(name) %>%
  summarise(
    n_hits       = sum(!is.na(p_match)),
    n_runs       = n_distinct(sample),
    mean_match   = mean(p_match, na.rm = TRUE),
    median_match = median(p_match, na.rm = TRUE),
    sd_match     = sd(p_match, na.rm = TRUE),
    se_match     = sd_match / sqrt(n_hits),
    .groups      = "drop"
  )

field_result <- read.csv("field_sampled_results.csv")

field_result

field_result_pct <- field_result %>%
  mutate(p_match = f_match * 100) %>%
  filter(sites != "Control")

field_clean <- field_result %>%
  mutate(p_match = as.numeric(f_match),
         intersect_bp = as.numeric(intersect_bp),
         unique_intersect_bp = as.numeric(unique_intersect_bp),
         ratio = unique_intersect_bp / intersect_bp
  ) %>%
  group_by(sample) %>%
  mutate(
    score = unique_intersect_bp / sum(unique_intersect_bp, na.rm = TRUE),
    like  = score * ratio
  ) %>%
  mutate(p_match = f_match * 100)

# filter score
field_clean <- field_clean %>%
  filter(unique_intersect_bp > 4000, score > 0.5, sites != "Control")

field_summary_stat <- field_clean %>%
  group_by(name)%>%
  summarise(
    n_hits       = sum(!is.na(like)),
    n_runs       = n_distinct(sample),
    mean_score   = mean(p_match, na.rm = TRUE),
    median_score = median(like, na.rm = TRUE),
    sd_score     = sd(like, na.rm = TRUE),
    se_score    = sd_score / sqrt(n_hits),
    .groups      = "drop"
  )

field_summary_stat

per_sample <- field_clean %>%
  group_by(name, sites) %>%
  summarise(like = mean(like, na.rm = TRUE), .groups = "drop")

# Global difference among species (non-parametric ANOVA)
kw <- kruskal.test( like ~ name, data = per_sample)
kw
