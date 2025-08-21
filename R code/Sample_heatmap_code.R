# Mapping sourmash matches for field samples
# ratio and score comparison of match confidance

rm(list=ls())

setwd("~/MRes Biosystematics/Project 2/data/results")
getwd()

library(tidyverse)
library(sf)
library(leaflet)
library(htmlwidgets)
library(viridisLite)
library(ggplot2)
library(dbplyr)
library(readr)

# load results data
field_result <- read.csv("Field_sampled_results.csv")


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


# filter score
field_clean <- field_clean %>%
  filter(unique_intersect_bp > 4000, score > 0.5, sites != "Control")

# compare between species 
con_summary_tbl <- field_clean %>%
  group_by(sites) %>%
  summarise(
    n_hits       = sum(!is.na(like)),                
    n_runs       = n_distinct(sites),                  
    mean_score   = mean(like, na.rm = TRUE),
    median_score = median(like, na.rm = TRUE),
    sd_score     = sd(like, na.rm = TRUE),           
    .groups = "drop"
  ) %>%
  mutate(
    se_score = sd_score / sqrt(n_hits))

#build heatmap
heatmap_df <- field_clean %>%
  filter(!is.na(like)) %>%
  select(sample, name, like) %>%
  distinct()

library(reshape2)

like_matrix <- reshape2::acast(heatmap_df, name ~ sample, value.var = "like", fill = 0)

heatmap(like_matrix, col = terrain.colors(256))

# Base R heatmap
heatmap(like_matrix, 
        Rowv = NA, Colv = NA,       
        scale = "none",             
        col = viridis::viridis(100), 
        margins = c(10, 10))

# pheatmap
library(pheatmap)

site_labels <- tibble(
  sample = c("barcode57_01", "barcode57_02", "barcode57_03",
             "barcode57_04", "barcode57_05", "barcode58_01", 
             "barcode59_01", "barcode59_02", "barcode59_03",
             "barcode59_04", "barcode59_05", "barcode60_01"),  
  site_name = c("Stock pond 1", "Stock pond 2", "Stock pond 3",
                "Stock pond 4", "Stock pond 5", "Bird Sancturary",
                "Men's pond 1", "Men's pond 2", "Men's pond 3",
                "Men's pond 4", "Men's pond 5", "Serpentine"))

heatmap_df_named <- heatmap_df %>%
  left_join(site_labels, by = "sample") %>%
  mutate(site_name = ifelse(is.na(site_name), sample, site_name))

# Convert to wide matrix
like_matrix <- heatmap_df_named %>%
  select(site_name, name, like) %>%
  pivot_wider(names_from = name, values_from = like, values_fill = 0) %>%
  column_to_rownames("site_name") %>%
  as.matrix()

like_matrix_transposed <- t(like_matrix)

pheatmap(
  like_matrix_transposed,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "white",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45  # Rotate column labels
)

# map binary observation vs sourmash
# Procambarus clarkii

binary_data <- read.csv("binary_procambarus_clarkii.csv")

binary_comparison <- binary_data %>%
  pivot_longer(cols = c(expected, sourmash), names_to = "type", values_to = "presence")

ggplot(binary_comparison, aes(x = type, y = sites, fill = factor(presence))) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                    name = "Estimate presence",
                    labels = c("Absent", "Present")) +
  coord_fixed(ratio = 1)+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

binary_data_match <- binary_data %>%
  mutate(match = case_when(
    expected == 1 & sourmash == 1 ~ "Expected & Present",
    expected == 1 & sourmash == 0 ~ "Expected Only",
    expected == 0 & sourmash == 1 ~ "Unexpected",
    TRUE ~ "Absent"
  ))

ggplot(binary_data_match, aes(x = species, y = sites, fill = match)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c(
    "Expected & Present" = "black",
    "Expected Only" = "grey40",
    "Unexpected" = "red",
    "Absent" = "white"
  )) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
