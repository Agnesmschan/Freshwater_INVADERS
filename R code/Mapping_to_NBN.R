# Mapping sourmash matches agaisnt NBN atlas and compare between species

rm(list=ls())

setwd("~/MRes Biosystematics/Project 2/data/results")
getwd()

library(tidyverse)
library(sf)
library(ggplot2)
library(viridisLite)
library(rnaturalearth)
library(rnaturalearthdata)
library(galah)
library(stringr)
library(scales)
library(ggrepel)


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

# select top match only
use_rank0_only <- TRUE

hits <- sourmash_result %>%
  mutate(f_match = suppressWarnings(as.numeric(f_match))) %>%
  { if (use_rank0_only) dplyr::filter(., gather_result_rank == 0) else . } %>%
  inner_join(meta_clean, by = c("sample" = "Run")) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(
    species = name,
    pmatch  = f_match * 100
  )

# ========= impatiens glandulifera =========
NBN_FILE <- "C:/Users/Agnes Chan/Documents/MRes Biosystematics/Project 2/data/results/nbn_impatiens_glandulifera.csv"
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
target_species <- "Impatiens glandulifera"
my_pts <- hits %>% filter(species == target_species)

if (nrow(my_pts) == 0) {
  stop("No sample hits found for '", target_species, "'. Check the 'name' column values.")
}

# ========= Basemap & extent =========
uk <- rnaturalearth::ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")

# ======== Plot Impatiens glandulifera top hits map ==========
Impatiens_glandulifera_map <- ggplot() +
  geom_sf(data = uk, fill = "grey98", colour = "grey80", linewidth = 0.5) +
  # NBN background points (faint)
  geom_point(data = nbn_sf, aes(x = longitude, y = latitude), size = 0.20, alpha = 0.10, colour = "coral4") +
  # Your samples with %match gradient
  geom_point(
    data = my_pts,
    aes(x = longitude, y = latitude, colour = pmatch),
    size = 2.8, alpha = 0.90
  ) +
  scale_colour_viridis_c("% match") +
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

Impatiens_glandulifera_map 

