library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(Rtsne)
library(umap)

# Set the correct file path
data <- read.csv("C:/Users/druschel/Downloads/D2vsD4_all genes 17000 for PCA.csv")

# Select only the expression columns for PCA (D2 and D4 samples)
expression_data <- data %>%
  dplyr::select(starts_with("D2_") | starts_with("D4_"))

# Transpose the data so that samples are rows and genes are columns
transposed_data <- as.data.frame(t(expression_data))

# Extract sample names and conditions from row names
sample_names <- rownames(transposed_data)
conditions <- ifelse(grepl("^D2_", sample_names), "D2", "D4")

# Log-transform the data (optional for RNA-seq data)
log_data <- log2(transposed_data + 1)

#### PCA - All ####
# Perform PCA
pca_result <- prcomp(log_data, center = TRUE, scale. = TRUE)

# Calculate variance explained
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var) * 100

# Prepare data for PCA plotting
pca_data <- as.data.frame(pca_result$x)
pca_data$Condition <- conditions
pca_data$Sample <- sample_names

# Plot PCA with sample names
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  labs(title = "PCA of All Samples (D2 and D4)",
       x = paste0("PC1: ", round(pca_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(pca_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### PCA - D2 Only ####
d2_data <- expression_data %>% dplyr::select(starts_with("D2_"))
d2_log_data <- log2(t(d2_data) + 1)
d2_pca_result <- prcomp(d2_log_data, center = TRUE, scale. = TRUE)
d2_var_explained <- d2_pca_result$sdev^2 / sum(d2_pca_result$sdev^2) * 100
d2_pca_data <- as.data.frame(d2_pca_result$x)
d2_pca_data$Sample <- colnames(d2_data)

# Plot D2 PCA with sample names
ggplot(d2_pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(size = 3) +
  labs(title = "PCA of D2 Samples Only",
       x = paste0("PC1: ", round(d2_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(d2_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### PCA - D4 Only ####
d4_data <- expression_data %>% dplyr::select(starts_with("D4_"))
d4_log_data <- log2(t(d4_data) + 1)
d4_pca_result <- prcomp(d4_log_data, center = TRUE, scale. = TRUE)
d4_var_explained <- d4_pca_result$sdev^2 / sum(d4_pca_result$sdev^2) * 100
d4_pca_data <- as.data.frame(d4_pca_result$x)
d4_pca_data$Sample <- colnames(d4_data)

# Plot D4 PCA with sample names
ggplot(d4_pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3, color = "red") +
  geom_text_repel(size = 3) +
  labs(title = "PCA of D4 Samples Only",
       x = paste0("PC1: ", round(d4_var_explained[1], 2), "% variance"),
       y = paste0("PC2: ", round(d4_var_explained[2], 2), "% variance")) +
  theme_minimal()

#### UMAP - All ####
umap_result <- umap(log_data, n_neighbors = 5, min_dist = 0.1)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Condition <- conditions
umap_df$Sample <- sample_names

# Plot UMAP with sample names
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  labs(title = "UMAP of All Samples") +
  theme_minimal()

#### UMAP - D2 Only ####
d2_umap_result <- umap(d2_log_data, n_neighbors = 2, min_dist = 0.1)
d2_umap_df <- as.data.frame(d2_umap_result$layout)
colnames(d2_umap_df) <- c("UMAP1", "UMAP2")
d2_umap_df$Sample <- colnames(d2_data)

# Plot UMAP for D2 samples with sample names
ggplot(d2_umap_df, aes(x = UMAP1, y = UMAP2, label = Sample)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(size = 3) +
  labs(title = "UMAP of D2 Samples Only") +
  theme_minimal()

#### UMAP - D4 Only ####
d4_umap_result <- umap(d4_log_data, n_neighbors = 2, min_dist = 0.1)
d4_umap_df <- as.data.frame(d4_umap_result$layout)
colnames(d4_umap_df) <- c("UMAP1", "UMAP2")
d4_umap_df$Sample <- colnames(d4_data)

# Plot UMAP for D4 samples with sample names
ggplot(d4_umap_df, aes(x = UMAP1, y = UMAP2, label = Sample)) +
  geom_point(size = 3, color = "red") +
  geom_text_repel(size = 3) +
  labs(title = "UMAP of D4 Samples Only") +
  theme_minimal()

#### t-SNE - All ####
tsne_result_all <- Rtsne(as.matrix(log_data), dims = 2, perplexity = 2)
tsne_df_all <- as.data.frame(tsne_result_all$Y)
colnames(tsne_df_all) <- c("tSNE1", "tSNE2")
tsne_df_all$Condition <- conditions
tsne_df_all$Sample <- sample_names

# Plot t-SNE with sample names
ggplot(tsne_df_all, aes(x = tSNE1, y = tSNE2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  labs(title = "t-SNE of All Samples") +
  theme_minimal()

#### t-SNE - D2 Only ####
d2_tsne_result <- Rtsne(as.matrix(d2_log_data), dims = 2, perplexity = 1)
d2_tsne_df <- as.data.frame(d2_tsne_result$Y)
colnames(d2_tsne_df) <- c("tSNE1", "tSNE2")
d2_tsne_df$Sample <- colnames(d2_data)

# Plot t-SNE for D2 samples with sample names
ggplot(d2_tsne_df, aes(x = tSNE1, y = tSNE2, label = Sample)) +
  geom_point(size = 3, color = "blue") +
  geom_text_repel(size = 3) +
  labs(title = "t-SNE of D2 Samples Only") +
  theme_minimal()

#### t-SNE - D4 Only ####
d4_tsne_result <- Rtsne(as.matrix(d4_log_data), dims = 2, perplexity = 1)
d4_tsne_df <- as.data.frame(d4_tsne_result$Y)
colnames(d4_tsne_df) <- c("tSNE1", "tSNE2")
d4_tsne_df$Sample <- colnames(d4_data)

# Plot t-SNE for D4 samples with sample names
ggplot(d4_tsne_df, aes(x = tSNE1, y = tSNE2, label = Sample)) +
  geom_point(size = 3, color = "red") +
  geom_text_repel(size = 3) +
  labs(title = "t-SNE of D4 Samples Only") +
  theme_minimal()
