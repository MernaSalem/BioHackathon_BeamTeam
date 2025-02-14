# Install and load required packages
if (!require("SNFtool")) install.packages("SNFtool")
if (!require("Rtsne")) install.packages("Rtsne")
if (!require("readr")) install.packages("readr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("cluster")) install.packages("cluster")

library(SNFtool)
library(Rtsne)
library(readr)
library(ggplot2)
library(pheatmap)
library(cluster)

# ---------------------- #
# Step 1: Load Data 
# ---------------------- #
methylation_data <- read_csv("methaylation.csv", col_names = TRUE)
expression_data <- read_csv("transcriptomics.csv", col_names = TRUE)

# Extract Feature Names
feature_names_meth <- methylation_data[[1]]  # CpG sites
feature_names_expr <- expression_data[[1]]   # Gene names

# Extract Sample Names from Column Headers
sample_names <- colnames(methylation_data)[-1]       
sample_names_expr <- colnames(expression_data)[-1]    

# Remove Feature Names Column & Convert to Matrix
methylation_matrix <- as.matrix(methylation_data[, -1])
expression_matrix <- as.matrix(expression_data[, -1])

# Assign Correct Row and Column Names BEFORE Transposing
rownames(methylation_matrix) <- feature_names_meth
rownames(expression_matrix) <- feature_names_expr
colnames(methylation_matrix) <- sample_names
colnames(expression_matrix) <- sample_names_expr

# Transpose the Data
T_Methylation_data <- t(methylation_matrix)
T_expression_data <- t(expression_matrix)

# Assign Sample Names as Row Names
rownames(T_Methylation_data) <- sample_names
rownames(T_expression_data) <- sample_names_expr

# ---------------------- #
# Step 2: Preprocessing 
# ---------------------- #
# Convert to numeric WITHOUT losing row names
mode(T_Methylation_data) <- "numeric"
mode(T_expression_data) <- "numeric"

# Get common samples
common_samples <- intersect(rownames(T_Methylation_data), rownames(T_expression_data))

# Subset data (ONLY if common_samples is not empty)
if (length(common_samples) > 0) {
  T_Methylation_data <- T_Methylation_data[common_samples, ]
  T_expression_data <- T_expression_data[common_samples, ]
} else {
  stop("No matching sample names found between datasets! Check row names before subsetting.")
}

# Verify that both matrices have the same samples
print(dim(T_Methylation_data))
print(dim(T_expression_data))

# Check for NA values and impute if necessary
imputeNA <- function(mat) {
  for (j in 1:ncol(mat)) {
    if (any(is.na(mat[, j]))) {
      mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
    }
  }
  return(mat)
}
T_Methylation_data <- imputeNA(T_Methylation_data)
T_expression_data <- imputeNA(T_expression_data)

# ---------------------- #
# Step 3: Compute Similarity Matrices 
# ---------------------- #
num_samples <- nrow(T_expression_data)
K_value <- min(10, num_samples - 1)

# Use original similarity calculation (Euclidean distance)
W_expr <- affinityMatrix(dist2(T_expression_data, T_expression_data), K = K_value, sigma = 0.5)
W_meth <- affinityMatrix(dist2(T_Methylation_data, T_Methylation_data), K = K_value, sigma = 0.5)

# Fuse Similarity Networks Using SNF
W_fused <- SNF(list(W_expr, W_meth), K = K_value, t = 20)

# Perform Spectral Clustering
num_clusters <- 3 
clusters <- spectralClustering(W_fused, num_clusters)

# ---------------------- #
# Step 4: Visualizations 
# ---------------------- #

# Hierarchical Clustering (with patient labels)
plot(hclust(as.dist(1 - W_fused)), labels = rownames(W_fused), main = "Hierarchical Clustering of SNF")

# PCA Plot with Patient Labels
pca_result <- prcomp(W_fused, scale = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Cluster = as.factor(clusters), Patient = rownames(W_fused))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Patient)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +  
  labs(title = "PCA of SNF Clustering", x = "PC1", y = "PC2") +
  theme_minimal()

# Heatmap of Fused Similarity Matrix with Cluster Labels
patient_clusters <- data.frame(Cluster = as.factor(clusters))
rownames(patient_clusters) <- rownames(W_fused)

pheatmap(W_fused, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_row = patient_clusters, 
         annotation_col = patient_clusters, 
         main = "Fused Similarity Matrix with Clusters")

# t-SNE Plot with Patient Labels
perplexity_value <- min(10, floor((num_samples - 1) / 3))  
if (perplexity_value > 1) {
  tsne_result <- Rtsne(W_fused, perplexity = perplexity_value)
  plot(tsne_result$Y, col = clusters, pch = 16, main = "t-SNE of SNF Clustered Patients")
  text(tsne_result$Y, labels = rownames(W_fused), pos = 4, cex = 0.7)  # Add patient labels
} else {
  print("Not enough samples for t-SNE, skipping this step.")
}

# ✅ 5️⃣ Silhouette Plot with Patient Labels
sil <- silhouette(clusters, dist(1 - W_fused))
plot(sil, main = "Silhouette Plot for SNF Clusters", border = NA)
text(sil[, 'sil_width'], labels = rownames(W_fused), pos = 3, cex = 0.7)  # Add patient labels
