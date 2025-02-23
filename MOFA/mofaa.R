
# Load Libraries

library(MOFA2)
library(tidyverse)
library(ggplot2)
library(readxl)
library(randomForest)
library(survminer)
library(tidyr)
library(dplyr)


#  Load Data

epigenomic_data <- read.csv("epigenomic_data_mat.csv", row.names = 1)
transcriptomic_data <- read.csv("Ntranscriptomic_data_mat.csv", row.names = 1)
metadata <- read_excel("Metad.xlsx")


# Ensure Column Names Match for Integration
# 
common_samples <- intersect(colnames(transcriptomic_data), colnames(epigenomic_data))
mRNA_counts <- transcriptomic_data[, common_samples]
Methylation_counts <- epigenomic_data[, common_samples]


# Scale Data

mRNA_scaled <- t(scale(t(log1p(mRNA_counts))))
Methylation_scaled <- t(scale(t(Methylation_counts), center = TRUE, scale = FALSE))

# Filter Low-Variance Genes

gene_variances <- apply(mRNA_scaled, 1, var, na.rm = TRUE)
gene_variances <- gene_variances[!is.na(gene_variances)]  # Remove NAs
high_var_genes <- gene_variances > quantile(gene_variances, 0.2, na.rm = TRUE)
mRNA_filtered <- mRNA_scaled[high_var_genes, ]

# Prepare MOFA Input

data_list <- list(mRNA = as.matrix(mRNA_filtered), Methylation = as.matrix(Methylation_scaled))
MOFAobject <- create_mofa(data_list)


# Set Model Options
# ========================
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
if ("sparsity" %in% names(model_opts)) model_opts$sparsity <- TRUE

# Training Options
train_opts <- get_default_training_options(MOFAobject)
MOFAobject_prepared <- prepare_mofa(MOFAobject, get_default_data_options(MOFAobject), model_opts, train_opts)

#  Run MOFA

MOFAobject_run <- run_mofa(MOFAobject_prepared, use_basilisk = TRUE)


# Visualization

plot_variance_explained(MOFAobject_run, plot_total = TRUE)
plot_variance_explained(MOFAobject_run, plot_total = FALSE)
plot_data_overview(MOFAobject_run)
plot_factor_cor(MOFAobject_run)

# ========================
# 🔗 Integrate Metadata
# ========================
metadata$group <- "group1"
MOFAobject_run@samples_metadata <- metadata

# Debug: Check Sample Name Mismatch
cat("Metadata Sample Names: ", metadata$sample, "\n")
cat("MOFA Sample Names: ", colnames(MOFAobject_run@samples_metadata), "\n")  #
# ensures that the sample names in your metadata and the MOFA object match perfectly for proper integration and analysis.

if (!all(metadata$sample %in% colnames(MOFAobject_run@samples_metadata))) {
  warning("Mismatch between metadata and MOFA sample names. Check and adjust metadata.")
}


# Correlate Factors with Covariates

correlate_factors_with_covariates(MOFAobject_run, 
                                  covariates = c("Sample_id_cpg",
                                                 "Sample_id_mrna",
                                                 "Sample_title_cpg",
                                                 "Sample_source_name_ch1_cpg",
                                                 "age",
                                                 "gender",
                                                 "sample"),
                                  plot = "log_pval")

#️ Clean and Map Factors

# Remove 'X' prefix from df_factors_wide
df_factors_wide <- get_factors(MOFAobject_run, as.data.frame = TRUE) %>%
  mutate(sample_clean = str_remove(sample, "^X"))

# Create a corrected sample mapping
sample_mapping <- data.frame(
  df_factors_sample = paste0("Sample_", 1:18),
  metadata_sample = c("1_L", "1_N", "2_L", "2_N", "3_L", "3_N", 
                      "4_L", "4_N", "5_L", "5_N", "6_L", "6_N", 
                      "7_L", "7_N", "9_L", "9_N", "10_L", "10_N")
) ####Build a mapping table that links the MOFA sample names (df_factors_sample) to  actual metadata sample names.
##
# Join factors with sample mapping
df_factors_mapped <- left_join(df_factors_wide, sample_mapping, 
                               by = c("sample_clean" = "metadata_sample"))#3Connect the MOFA factor data (df_factors_wide) with sample mapping table.

# Join mapped factors with metadata
df_plot <- left_join(df_factors_mapped, metadata, by = c("sample_clean" = "sample")) #Add the metadata information (like age, gender, etc.) to the mapped factors.
##now w have  mofa factors+sample mapping+metadata

# Fill NAs in Sample_source_name_ch1_cpg with 'Unknown'
df_plot$Sample_source_name_ch1_cpg[is.na(df_plot$Sample_source_name_ch1_cpg)] <- "Unknown"

# Pivot to Wide Format for Plotting

df_plot_wide <- df_plot %>%
  pivot_wider(names_from = factor, values_from = value) ##for spreading the factor values across new columns.

# Check new column names
print(colnames(df_plot_wide))


#  Plot Factor1 vs Factor3

ggplot(df_plot_wide, aes(x = Factor1, y = Factor3, color = Sample_source_name_ch1_cpg)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("MOFA: Factor 1 vs Factor 3") +
  labs(color = "Sample Source") +
  theme(plot.title = element_text(hjust = 0.5))


#  Extract Top Features for Factors

extract_top_features <- function(factor_num, view) {
  weights <- get_weights(MOFAobject_run, view = view, factors = factor_num, as.data.frame = TRUE)
  weights %>%
    filter(!is.na(value)) %>%
    arrange(desc(abs(value))) %>%
    slice(1:10)
}

# Extract top features
top_features_f1_mRNA <- extract_top_features(1, "mRNA")
top_features_f1_methyl <- extract_top_features(1, "Methylation")
top_features_f6_mRNA <- extract_top_features(6, "mRNA")
top_features_f6_methyl <- extract_top_features(6, "Methylation")


# Plot Top Features

plot_top_features <- function(data, title) {
  ggplot(data, aes(x = reorder(feature, value), y = value, fill = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = title, x = "Features", y = "Weight") +
    theme(plot.title = element_text(hjust = 0.5))
}

# Plot top features
plot_top_features(top_features_f1_mRNA, "Top 10 Features for Factor 1 (mRNA)")
plot_top_features(top_features_f1_methyl, "Top 10 Features for Factor 1 (Methylation)")
plot_top_features(top_features_f6_mRNA, "Top 10 Features for Factor 6 (mRNA)")
plot_top_features(top_features_f6_methyl, "Top 10 Features for Factor 6 (Methylation)")
