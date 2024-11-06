# Load necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(umap)
library(Rtsne)

#### User Inputs ####

# Define the number of groups (e.g., 2 for D2 and D4)
num_groups <- 2

# Define the names of each group
group_names <- c("D2", "D4")  # Adjust or input as needed

#### Data Preparation ####

# Read in the data
data <- read.csv("D2vsD4_all genes 17000 for PCA.csv")

# Select only expression columns based on group names
expression_data <- data %>%
  dplyr::select(matches(paste(group_names, collapse = "|")))

# Transpose the data so that samples are rows and genes are columns
transposed_data <- as.data.frame(t(expression_data))

# Log-transform the data (optional for RNA-seq data)
log_data <- log2(transposed_data + 1)

#### PCA - All Samples Combined ####

# Perform PCA
pca_result <- prcomp(log_data, center = TRUE, scale. = TRUE)

# Check variance explained by each principal component
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var) * 100

# Prepare data for plotting
pca_data <- as.data.frame(pca_result$x)

# Assign group labels dynamically based on sample counts in each group
sample_counts <- sapply(group_names, function(g) sum(grepl(g, colnames(expression_data))))
pca_data$Condition <- unlist(mapply(rep, group_names, sample_counts))

# Plot PCA for all samples
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of All Samples",
       x = paste0("PC1: ", round(pca_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### PCA for Each Group Separately ####

# Set a small variance threshold to filter out low-variance columns
variance_threshold <- 1e-6

for (group in group_names) {
  # Subset and transpose data for the current group
  group_data <- expression_data %>% dplyr::select(starts_with(group))
  group_transposed <- as.data.frame(t(group_data))
  group_log_data <- log2(group_transposed + 1)
  
  # Filter out columns with near-zero variance
  group_log_data <- group_log_data[, apply(group_log_data, 2, var) > variance_threshold]
  
  # Perform PCA on the current group
  group_pca_result <- prcomp(group_log_data, center = TRUE, scale. = TRUE)
  group_pca_data <- as.data.frame(group_pca_result$x)
  group_var_explained <- group_pca_result$sdev^2 / sum(group_pca_result$sdev^2) * 100
  
  # Plot PCA for the current group
  ggplot(group_pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = ifelse(group == group_names[1], "blue", "red")) +
    labs(title = paste("PCA of", group, "Samples Only"),
         x = paste0("PC1: ", round(group_var_explained[1], 2), "% variance"),
         y = paste0("PC2: ", round(group_var_explained[2], 2), "% variance")) +
    theme_minimal()
}

#### UMAP - All Samples Combined ####

# Determine n_neighbors for UMAP on all samples
total_samples <- sum(sample_counts)
umap_neighbors <- max(1, total_samples - 1)  # Ensure n_neighbors is always less than sample count
umap_result <- umap(log_data, n_neighbors = umap_neighbors, min_dist = 0.1)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Condition <- unlist(mapply(rep, group_names, sample_counts))

# Plot UMAP for all samples
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "UMAP of All Samples") +
  theme_minimal()

#### UMAP for Each Group Separately ####

for (group in group_names) {
  # Subset and transpose data for the current group
  group_data <- expression_data %>% dplyr::select(starts_with(group))
  group_transposed <- as.data.frame(t(group_data))
  group_log_data <- log2(group_transposed + 1)
  
  # Calculate n_neighbors based on the number of samples in the current group (rows)
  group_samples <- nrow(group_log_data)
  group_neighbors <- max(1, group_samples - 1)  # Ensures n_neighbors is always less than sample size
  
  # Run UMAP for the current group
  group_umap_result <- umap(group_log_data, n_neighbors = group_neighbors, min_dist = 0.1)
  group_umap_df <- as.data.frame(group_umap_result$layout)
  colnames(group_umap_df) <- c("UMAP1", "UMAP2")
  
  # Plot and print UMAP for the current group
  umap_plot <- ggplot(group_umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 3, color = ifelse(group == group_names[1], "blue", "red")) +
    labs(title = paste("UMAP of", group, "Samples Only")) +
    theme_minimal()
  print(umap_plot)
}

#### t-SNE - All Samples Combined ####

# Determine perplexity for t-SNE (suggested: between 1/3 and half of total samples)
tsne_perplexity <- 2

# Run t-SNE on all samples with determined perplexity
tsne_result_all <- Rtsne(as.matrix(log_data), dims = 2, perplexity = tsne_perplexity)
tsne_df_all <- as.data.frame(tsne_result_all$Y)
colnames(tsne_df_all) <- c("tSNE1", "tSNE2")
tsne_df_all$Condition <- unlist(mapply(rep, group_names, sample_counts))

# Plot t-SNE for all samples
ggplot(tsne_df_all, aes(x = tSNE1, y = tSNE2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of All Samples") +
  theme_minimal()

#### t-SNE for Each Group Separately ####

#### t-SNE for Each Group Separately ####

for (group in group_names) {
  # Subset and transpose data for the current group
  group_data <- expression_data %>% dplyr::select(starts_with(group))
  group_transposed <- as.data.frame(t(group_data))
  group_log_data <- log2(group_transposed + 1)
  
  # Calculate perplexity based on sample size for current group (rows)
  group_samples <- nrow(group_log_data)
  group_perplexity <- 1
    # didn;t work, try adding back in: max(1, group_samples - 1)  # Ensures perplexity is always less than sample size
  
  # Run t-SNE for the current group
  group_tsne_result <- Rtsne(as.matrix(group_log_data), dims = 2, perplexity = group_perplexity)
  group_tsne_df <- as.data.frame(group_tsne_result$Y)
  colnames(group_tsne_df) <- c("tSNE1", "tSNE2")
  
  # Plot and print t-SNE for the current group
  tsne_plot <- ggplot(group_tsne_df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(size = 3, color = ifelse(group == group_names[1], "blue", "red")) +
    labs(title = paste("t-SNE of", group, "Samples Only")) +
    theme_minimal()
  print(tsne_plot)
}

