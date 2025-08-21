# Mapping sourmash matches and compare between groups
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
  ungroup()

# filter score
results_clean <- results_clean %>%
  filter(unique_intersect_bp > 4000, score > 0.5)

# merge tables by accessions
sourmash_hits <- results_clean%>%
  inner_join(meta_clean, by = c("sample" = "Run")) %>%
  filter(!is.na(latitude), !is.na(longitude))

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
pal <- colorNumeric(viridis(256, direction = -1), domain = pts$like, na.color = "transparent")

popup = ~sprintf(
  "<b>%s</b><br/>Run: %s<br/>%% match: %.4f<br/>Score: %.3f<br/>Ratio: %.3f<br/>Like: %.3f<br/>Country: %s",
  name, sample, pmatch, score, ratio, like, as.character(geo_loc_name_country)
)

# Update the map function to use `like` for fillColor and legend
make_map <- function(dat_sf, title, outfile) {
  m <- leaflet(dat_sf) %>%
    addTiles() %>%
    addCircleMarkers(
      lng = ~longitude, lat = ~latitude,
      radius = 6, stroke = FALSE, fillOpacity = 0.85,
      fillColor = ~pal(like),  # <- use 'like' for coloring
      clusterOptions = markerClusterOptions(),
      popup = popup  # <- use updated popup
    ) %>%
    addLegend("bottomright", pal = pal, values = ~like,
              title = paste0(title, "<br/>Confidence index"),
              labFormat = labelFormat(), opacity = 1)
  htmlwidgets::saveWidget(m, outfile, selfcontained = TRUE)
  m
}

# Generate maps
map_plants_confidence <- make_map(pts %>% filter(group == "Plant"), "Plants", "map_plants.html")
map_invertebrates_confidence <- make_map(pts %>% filter(group == "Invertebrate"), "Invertebrates", "map_invertebrates.html")

# View map
map_plants_confidence
map_invertebrates_confidence

# compare between species 
con_summary_tbl <- hits_pi %>%
  group_by(name) %>%
  summarise(
    n_hits       = sum(!is.na(like)),                
    n_runs       = n_distinct(sample),                  
    mean_score   = mean(like, na.rm = TRUE),
    median_score = median(like, na.rm = TRUE),
    sd_score     = sd(like, na.rm = TRUE),           
    .groups = "drop"
  ) %>%
  mutate(
    se_score = sd_score / sqrt(n_hits))

# compare group between groups
group_con_summary_tbl <- hits_pi %>%
  group_by(group) %>%
  summarise(
    n_hits       = sum(!is.na(like)),                
    n_runs       = n_distinct(sample),                  
    mean_score   = mean(like, na.rm = TRUE),
    median_score = median(like, na.rm = TRUE),
    sd_score     = sd(like, na.rm = TRUE),           
    .groups = "drop"
  ) %>%
  mutate(
    se_score = sd_score / sqrt(n_hits))

# Boxplot + jitter of %match by group
p_box <- ggplot(hits_pi, aes(x = name, y = like)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), name = "Confidence") +
  xlab(NULL) +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

p_box

# Non-parametric test (Wilcoxon) and effect size (Cliff's delta)
comp_dat <- hits_pi %>% filter(group %in% c("Plant","Invertebrate"))
wilc <- wilcox.test(like ~ group, data = comp_dat, exact = FALSE)
cdel <- effsize::cliff.delta(like ~ group, data = comp_dat)

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
  filter(name %in% species_to_test, is.finite(like))

# collapse to one observation per sample×species (use mean() if you prefer)
per_sample <- species_comp_dat %>%
  group_by(name, sample) %>%
  summarise(like = max(like, na.rm = TRUE), .groups = "drop")

# Global difference among species (non-parametric ANOVA)
kw <- kruskal.test(like ~ name, data = per_sample)
kw

# Pairwise Wilcoxon (BH-adjusted p-values)
pw <- pairwise.wilcox.test(per_sample$like, per_sample$name,
                           p.adjust.method = "BH", exact = FALSE)
pw

#build heatmap
heatmap_df <- sourmash_hits %>%
  filter(!is.na(like)) %>%
  select(sample, name, like) %>%
  distinct()
# Add an index to help filter every nth label
heatmap_df <- heatmap_df %>%
  mutate(sample_index = as.integer(factor(sample)))

# Choose how often to label (e.g. every 10th)
every_n <- 10
sample_labels <- heatmap_df %>%
  distinct(sample_index, sample) %>%
  arrange(sample_index) %>%
  filter(sample_index %% every_n == 1)

SRA_heatmap <- ggplot(heatmap_df, aes(x = sample_index, y = name, fill = like)) +
  geom_tile(color = "white") +
  scale_x_continuous(
    breaks = sample_labels$sample_index,
    labels = sample_labels$sample,
    limits = c(0, NA), expand = c(0, 0)
  ) +
  scale_fill_viridis_c(direction = -1, name = "Confidence Score", na.value = "grey90") +
  labs(x = "Subset Sample", y = "Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

SRA_heatmap


# ========= static map with high confidence (per species) =========
# change map parameters for each species
NBN_FILE <- "C:/Users/Agnes Chan/Documents/MRes Biosystematics/Project 2/data/results/nbn_impatiens_glandulifera.csv" #change file for each 
file.exists(NBN_FILE)

nbn_raw <- readr::read_delim(NBN_FILE, delim = ifelse(grepl("\\.tsv$", NBN_FILE, ignore.case=TRUE), "\t", ","),
                             show_col_types = FALSE, guess_max = 100000)

# Try decimal lat/lon first (Darwin Core), else easting/northing (OSGB)
if (all(c("Latitude (WGS84)","Longitude (WGS84)") %in% names(nbn_raw))) {
  # Darwin Core-style CSV from NBN with WGS84 names
  nbn_tbl <- nbn_raw %>%
    mutate(
      latitude  = suppressWarnings(as.numeric(`Latitude (WGS84)`)),
      longitude = suppressWarnings(as.numeric(`Longitude (WGS84)`))
    ) %>%
    filter(is.finite(latitude), is.finite(longitude))
  
  nbn_sf <- sf::st_as_sf(nbn_tbl, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
  
} else if (all(c("decimalLatitude","decimalLongitude") %in% names(nbn_raw))) {
  # Alternative Darwin Core naming
  nbn_tbl <- nbn_raw %>%
    rename(latitude = decimalLatitude, longitude = decimalLongitude) %>%
    mutate(across(c(latitude, longitude), ~ suppressWarnings(as.numeric(.)))) %>%
    filter(is.finite(latitude), is.finite(longitude))
  
  nbn_sf <- sf::st_as_sf(nbn_tbl, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
  
} else if (all(c("easting","northing") %in% names(nbn_raw))) {
  # British National Grid -> WGS84
  nbn_tbl <- nbn_raw %>%
    mutate(
      easting  = suppressWarnings(as.numeric(easting)),
      northing = suppressWarnings(as.numeric(northing))
    ) %>%
    filter(is.finite(easting), is.finite(northing))
  
  nbn_sf <- sf::st_as_sf(nbn_tbl, coords = c("easting","northing"), crs = 27700) %>%
    sf::st_transform(4326) %>%
    mutate(
      longitude = sf::st_coordinates(.)[,1],
      latitude  = sf::st_coordinates(.)[,2]
    )
  
} else {
  stop("Couldn't find coordinates. Expected {Latitude (WGS84), Longitude (WGS84)} or {decimalLatitude, decimalLongitude} or {easting, northing}.")
}

# ========= Filter your samples to the species of interest =========
target_species <- "Impatiens glandulifera" #change species for each map
my_pts <- sourmash_hits %>% filter(name == target_species)

if (nrow(my_pts) == 0) {
  stop("No sample hits found for '", target_species, "'. Check the 'name' column values.")
}

# ========= Basemap & extent =========
uk <- rnaturalearth::ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")

# ======== Plot Impatiens glandulifera top hits map ==========
Impatiens_glandulifera_map <- ggplot() +
  geom_sf(data = uk, fill = "grey98", colour = "grey80", linewidth = 0.5) +
  # NBN background points (faint)
  geom_point(data = nbn_sf, aes(x = longitude, y = latitude), size = 0.2, alpha = 0.1, colour = "coral4") +
  # Your samples with %match gradient
  geom_point(
    data = my_pts,
    aes(x = longitude, y = latitude, colour = like),
    size = 2.8, alpha = 0.7
  ) +
  scale_colour_viridis_c("Confidence index", direction = -1) +
  coord_sf(
    datum = NA,
    xlim = st_bbox(uk)[c("xmin","xmax")],
    ylim = st_bbox(uk)[c("ymin","ymax")]
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title = element_text(face = "bold")
  )

# Map 1 - species 1
Impatiens_glandulifera_map

# Map 2 - species 2
Procambarus_clarkii_map

# Map 3 - species 3
Pontederia_crassipes_map

# Map 4 - species 4
Dreissena_polymorpha_map
