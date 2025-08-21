# Mapping sourmash matches and compare between groups
# raw similarity comparison of f_match

rm(list=ls())

setwd("~/MRes Biosystematics/Project 2/data/results")
getwd()

# load packages
library(tidyverse)
library(sf)
library(leaflet)
library(htmlwidgets)
library(viridisLite)
library(ggplot2)
library(scales)
library(effsize)
library(taxize)

# load results data
sourmash_result <- read.csv("SRA_search_results.csv")
metadata <- read.csv("Freshwater_metadata.csv")

# Clean up table
meta_clean <- metadata %>%
  rename(
    latitude  = `geographic_location_.latitude.`,
    longitude = `geographic_location_.longitude.`
  ) %>%
  mutate(
    latitude  = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  )

results_clean <- sourmash_result %>%
  mutate(f_match = suppressWarnings(as.numeric(f_match)))

# merge tables by accessions
sourmash_hits <- sourmash_result %>%
  inner_join(meta_clean, by = c("sample" = "Run")) %>%
  filter(!is.na(latitude), !is.na(longitude))

library(readr)
library(dplyr)

map_tbl <- readr::read_csv("group_map.csv", show_col_types = FALSE) %>%
  mutate(
    group = case_when(
      str_to_lower(group) %in% c("plant","plants") ~ "Plant",
      str_to_lower(group) %in% c("invertebrate","invertebrates","invert") ~ "Invertebrate",
      TRUE ~ "Other"
    )
  ) %>%
  distinct(name, .keep_all = TRUE)

sourmash_hits <- sourmash_hits %>%
  left_join(map_tbl %>% select(name, group), by = "name")


hits_pi <- sourmash_hits %>% filter(group %in% c("Plant","Invertebrate")) %>%
  mutate(pmatch = f_match * 100)

pts <- st_as_sf(hits_pi, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
pal <- colorNumeric(viridis(256), domain = pts$pmatch, na.color = "transparent")


make_map <- function(dat_sf, title, outfile) {
  m <- leaflet(dat_sf) %>%
    addTiles() %>%
    addCircleMarkers(
      lng = ~longitude, lat = ~latitude,
      radius = 6, stroke = FALSE, fillOpacity = 0.85,
      fillColor = ~pal(pmatch),
      clusterOptions = markerClusterOptions(),
      popup = ~sprintf(
        "<b>%s</b><br/>Run: %s<br/>%% match: %.2f<br/>Rank: %s<br/>Country: %s",
        name, sample, pmatch, as.character(gather_result_rank), as.character(geo_loc_name_country)
      )
    ) %>%
    addLegend("bottomright", pal = pal, values = ~pmatch,
              title = paste0(title, "<br/>% match"),
              labFormat = labelFormat(suffix = "%"), opacity = 1)
  htmlwidgets::saveWidget(m, outfile, selfcontained = TRUE)
  m
}

map_plants        <- make_map(pts %>% dplyr::filter(group == "Plant"),       "Plants",        "map_plants.html")
map_invertebrates <- make_map(pts %>% dplyr::filter(group == "Invertebrate"),"Invertebrates", "map_invertebrates.html")

map_plants
map_invertebrates

# compare between species 
summary_tbl <- hits_pi %>%
  group_by(name) %>%
  summarise(
    n_hits       = sum(!is.na(f_match)),                
    n_runs       = n_distinct(sample),                  
    mean_match   = mean(f_match, na.rm = TRUE),
    median_match = median(f_match, na.rm = TRUE),
    sd_match     = sd(f_match, na.rm = TRUE),           
    .groups = "drop"
  ) %>%
  mutate(
    se_match = sd_match / sqrt(n_hits),                 
    # convert all to percent
    across(c(mean_match, median_match, sd_match, se_match), ~ .x * 100)
  )

# compare group between groups
group_summary_tbl <- hits_pi %>%
  group_by(group) %>%
  summarise(
    n_hits       = sum(!is.na(f_match)),                
    n_runs       = n_distinct(sample),                  
    mean_match   = mean(f_match, na.rm = TRUE),
    median_match = median(f_match, na.rm = TRUE),
    sd_match     = sd(f_match, na.rm = TRUE),           
    .groups = "drop"
  ) %>%
  mutate(
    se_match = sd_match / sqrt(n_hits),                 
    # convert all to percent
    across(c(mean_match, median_match, sd_match, se_match), ~ .x * 100)
  )

print(summary_tbl)
readr::write_csv(summary_tbl, "species_summary.csv")
readr::write_csv(group_summary_tbl, "group_summary.csv")

# Boxplot + jitter of %match by group
p_box <- ggplot(hits_pi, aes(x = name, y = f_match)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), name = "% match") +
  xlab(NULL) +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

p_box

# Non-parametric test (Wilcoxon) and effect size (Cliff's delta)
comp_dat <- hits_pi %>% filter(group %in% c("Plant","Invertebrate"))
wilc <- wilcox.test(f_match ~ group, data = comp_dat, exact = FALSE)
cdel <- effsize::cliff.delta(f_match ~ group, data = comp_dat)

cat("\nWilcoxon rank-sum test:\n")
print(wilc)

cat("\nCliff's delta (effect size):\n")
print(cdel)

# Non-parametric test by species
species_to_test <- c(
  "Azolla filiculoides","Corbicula fluminalis","Dreissena polymorpha",
  "Hottonia palustris","Impatiens glandulifera","Lemna minuta",
  "Limnoperna fortunei","Pontederia crassipes","Procambarus clarkii",
  "Sinanodonta woodiana"
)

species_comp_dat <- hits_pi %>%
  filter(name %in% species_to_test, is.finite(f_match))

# collapse to one observation per sample×species (use mean() if you prefer)
per_sample <- species_comp_dat %>%
  group_by(name, sample) %>%
  summarise(f_match = max(f_match, na.rm = TRUE), .groups = "drop")

# Global difference among species (non-parametric ANOVA)
kw <- kruskal.test(f_match ~ name, data = per_sample)
kw

# Pairwise Wilcoxon (BH-adjusted p-values)
pw <- pairwise.wilcox.test(per_sample$f_match, per_sample$name,
                           p.adjust.method = "BH", exact = FALSE)
pw


# 1) Summary stats per group (mean, SD, SE) across samples
summary_counts <- per_run %>%
  group_by(group) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_n = mean(n),                 # mean hits per sample
    sd_n   = sd(n),                   # SD across samples
    se_n   = sd_n / sqrt(n_samples),  # SE of the mean
    .groups = "drop"
  )

summary_counts

# 2) Paired Wilcoxon test on per-sample counts (Plant vs Invertebrate)
wide_counts <- per_run %>%
  select(sample, group, n) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

wilc_paired <- wilcox.test(wide_counts$Plant, wide_counts$Invertebrate,
                           paired = TRUE, exact = FALSE)

p_label <- if (wilc_paired$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", wilc_paired$p.value)

# 3) Bar plot with SE error bars + colours
ymax <- max(summary_counts$mean_n + summary_counts$se_n) * 1.15

p_counts <- ggplot(summary_counts, aes(x = group, y = mean_n, fill = group)) +
  geom_col(width = 0.6, colour = "black") +
  geom_errorbar(aes(ymin = mean_n - mean_n, ymax = mean_n + se_n),
                width = 0.15, linewidth = 0.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    x = NULL, y = "Mean hits per sample",
    subtitle = p_label
  ) +
  coord_cartesian(ylim = c(0, ymax)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

p_counts

paired_counts <- ggplot(per_run, aes(x = group, y = n, group = sample, colour = group)) +
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +
  geom_line(aes(group = sample), colour = "grey70", alpha = 0.3) +
  scale_colour_brewer(palette = "Paired") +
  labs(x = NULL, y = "Hits per sample", title = "Per-sample counts (paired)") +
  theme_minimal() + theme(legend.position = "none")

library(forcats)
library(multcompView)

