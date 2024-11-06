# Load necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(umap)
library(Rtsne)

#### Prep Data ####

# Read in the data
data <- read.csv("D2vsD4_all genes 17000 for PCA.csv")

# Select only the expression columns for PCA (D2 and D4 samples)
expression_data <- data %>%
  dplyr::select(starts_with("D2_") | starts_with("D4_"))

# Transpose the data so that samples are rows and genes are columns
transposed_data <- as.data.frame(t(expression_data))

# Log-transform the data (optional for RNA-seq data)
log_data <- log2(transposed_data + 1)

#### PCA - All ####
# Perform PCA
pca_result <- prcomp(log_data, center = TRUE, scale. = TRUE)

# Check variance explained by each principal component
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var) * 100
print(pca_var_explained)  # Optional: view variance explained

# Prepare data for plotting
pca_data <- as.data.frame(pca_result$x)
pca_data$Condition <- c(rep("D2", 4), rep("D4", 5))  # Labeling 4 samples for D2 and 5 for D4

# PCA for all samples (already in the existing script)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of All Samples (D2 and D4)",
       x = paste0("PC1: ", round(pca_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### PCA for each group ####

# Set a small variance threshold to filter out low-variance columns
variance_threshold <- 1e-6

# Subset the data for D2 samples
d2_data <- expression_data %>% dplyr::select(starts_with("D2_"))
# Filter out columns with near-zero variance in D2 data
d2_log_data <- d2_log_data[, apply(d2_log_data, 2, var) > variance_threshold]

# Perform PCA on D2 samples with filtered data
d2_pca_result <- prcomp(d2_log_data, center = TRUE, scale. = TRUE)

# Prepare data for D2 PCA plot
d2_pca_data <- as.data.frame(d2_pca_result$x)
d2_var_explained <- d2_pca_result$sdev^2 / sum(d2_pca_result$sdev^2) * 100

# Plot PCA for D2 samples only
ggplot(d2_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "PCA of D2 Samples Only",
       x = paste0("PC1: ", round(d2_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(d2_var_explained[2], 2), "% variance")) +
  theme_minimal()

# Subset the data for D4 samples
d4_data <- expression_data %>% dplyr::select(starts_with("D4_"))
# Filter out columns with near-zero variance in D4 data
d4_log_data <- d4_log_data[, apply(d4_log_data, 2, var) > variance_threshold]

# Perform PCA on D4 samples with filtered data
d4_pca_result <- prcomp(d4_log_data, center = TRUE, scale. = TRUE)

# Prepare data for D4 PCA plot
d4_pca_data <- as.data.frame(d4_pca_result$x)
d4_var_explained <- d4_pca_result$sdev^2 / sum(d4_pca_result$sdev^2) * 100

# Plot PCA for D4 samples only
ggplot(d4_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, color = "red") +
  labs(title = "PCA of D4 Samples Only",
       x = paste0("PC1: ", round(d4_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(d4_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### UMAP - All ####

# Prepare data for UMAP by using the scaled data used for PCA (without labels)
umap_data <- log_data  # Assuming 'log_data' contains the scaled, transformed values

# Run UMAP on the transformed data, used n_neighbors of 5 but check
umap_result <- umap(umap_data, n_neighbors = 5, min_dist = 0.1)

# Convert the UMAP output to a data frame
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Add the sample condition information for plotting
umap_df$Condition <- rep(c("D2", "D4"), times = c(4, 5))  # Adjust based on actual sample counts

# Plot UMAP
library(ggplot2)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "UMAP of Gene Expression Data") +
  theme_minimal()

#### UMAP - D2 and D4 ####
# Subset for D2 samples
d2_data <- expression_data %>% dplyr::select(starts_with("D2_"))
d2_transposed <- as.data.frame(t(d2_data))
d2_log_data <- log2(d2_transposed + 1)

# Run UMAP for D2 samples
d2_umap_result <- umap(d2_log_data, n_neighbors = 2, min_dist = 0.1)
d2_umap_df <- as.data.frame(d2_umap_result$layout)
colnames(d2_umap_df) <- c("UMAP1", "UMAP2")

# Plot UMAP for D2 samples
ggplot(d2_umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "UMAP of D2 Samples Only") +
  theme_minimal()

# Subset for D4 samples
d4_data <- expression_data %>% dplyr::select(starts_with("D4_"))
d4_transposed <- as.data.frame(t(d4_data))
d4_log_data <- log2(d4_transposed + 1)

# Run UMAP for D4 samples
d4_umap_result <- umap(d4_log_data, n_neighbors = 2, min_dist = 0.1)
d4_umap_df <- as.data.frame(d4_umap_result$layout)
colnames(d4_umap_df) <- c("UMAP1", "UMAP2")

# Plot UMAP for D4 samples
ggplot(d4_umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 3, color = "red") +
  labs(title = "UMAP of D4 Samples Only") +
  theme_minimal()

#### tSNE - All ####

# Run t-SNE on all samples with a lower perplexity
tsne_result_all <- Rtsne(as.matrix(tsne_data), dims = 2, perplexity = 2)  # Adjusted perplexity
tsne_df_all <- as.data.frame(tsne_result_all$Y)
colnames(tsne_df_all) <- c("tSNE1", "tSNE2")
tsne_df_all$Condition <- rep(c("D2", "D4"), times = c(4, 5))  # Adjust for sample counts

# Plot t-SNE for all samples
ggplot(tsne_df_all, aes(x = tSNE1, y = tSNE2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of All Samples") +
  theme_minimal()

#### tSNE - D2/D4 Separate ####

# Run t-SNE for D2 samples with a very low perplexity
d2_tsne_result <- Rtsne(as.matrix(d2_log_data), dims = 2, perplexity = 1)
d2_tsne_df <- as.data.frame(d2_tsne_result$Y)
colnames(d2_tsne_df) <- c("tSNE1", "tSNE2")

# Plot t-SNE for D2 samples
ggplot(d2_tsne_df, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "t-SNE of D2 Samples Only") +
  theme_minimal()

# Run t-SNE for D4 samples with a very low perplexity
d4_tsne_result <- Rtsne(as.matrix(d4_log_data), dims = 2, perplexity = 1)
d4_tsne_df <- as.data.frame(d4_tsne_result$Y)
colnames(d4_tsne_df) <- c("tSNE1", "tSNE2")

# Plot t-SNE for D4 samples
ggplot(d4_tsne_df, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 3, color = "red") +
  labs(title = "t-SNE of D4 Samples Only") +
  theme_minimal()