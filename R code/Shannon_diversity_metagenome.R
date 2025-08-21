# Shannon diversity index for metagenome

setwd("~/MRes Biosystematics/Project 2/data")
getwd()

library(vegan)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

metagenome_data <- read_tsv("marti_assignments_all_levels_2025_v2.tsv", show_col_types = FALSE)

print(colnames(metagenome_data))

count_data <- metagenome_data %>%
  select(starts_with("AGNES"))

count_matrix <- as.data.frame(t(count_data))

colnames(count_matrix) <- metagenome_data$Name

count_matrix[] <- lapply(count_matrix, as.numeric)

count_matrix[is.na(count_matrix)] <- 0

alpha_shannon <- diversity(count_matrix, index = "shannon")
print("Shannon Alpha Diversity:")
print(alpha_shannon)

beta_bray <- vegdist(count_matrix, method = "bray")
print("Bray-Curtis Beta Diversity Metrix:")
print(as.matrix(beta_bray))

pooled_counts <- colSums(count_matrix)
gamma_shannon <- diversity(pooled_counts, index = "shannon")
cat("Gamma Shannon Diversity:", gamma_shannon, "\n")

# PCoA for beta diversity
pcoa_result <- cmdscale(beta_bray, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(pcoa_df)

pcoa_df$Barcode <- gsub(" \\(.*", "", pcoa_df$Sample)

site_labels <- tibble(
  Barcode = c("AGNES057", "AGNES058", "AGNES059", "AGNES060", "AGNES061"),  
  site_name = c("Stock pond", "Bird Sancturary", "Men's pond", "Serpentine", "Control"))

pcoa_df <- pcoa_df %>%
  left_join(site_labels, by = "Barcode") %>%
  mutate(site_name = ifelse(is.na(site_name), Barcode, site_name))

eig_vals <- pcoa_result$eig
variance_explained <- eig_vals / sum(eig_vals)
cat("Variance explained:\n")
print(variance_explained[1:2]) 


beta_PCoA <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = site_name)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_text_repel(size = 3) +
  labs(x = paste0("PCoA1 (", round(variance_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(variance_explained[2] * 100, 1), "%)"),
       color = "Sites") +
  theme_minimal()

beta_PCoA

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = site_name, label = site_name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = Inf) +
  labs(x = paste0("PCoA1 (", round(variance_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(variance_explained[2] * 100, 1), "%)"),
       color = "Sites") +
  theme_minimal()
