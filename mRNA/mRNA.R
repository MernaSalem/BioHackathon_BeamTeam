# Load necessary libraries
library(DESeq2)
library(readxl)
library(pheatmap)
library(ggplot2)
library(readr)
library(NMF)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
# loading data 
pheno <- read.xlsx("pheno.xlsx", sheet = 1, rowNames = TRUE)
data_file <- read.xlsx("GSE224377_raw_counts_GRCh38.p13_NCBI.xlsx",sheet=1, rowNames = TRUE)
data_matrix <- as.matrix(data_file) 

# Remove rows where all values r  zero
filtered_data <- data_matrix[rowSums(data_matrix != 0) > 0, ]

# Save filtered data
write.csv(filtered_data,"GSE224377_filtered.csv", row.names=TRUE)
#load daata
Raw_counts <- read.csv("GSE224377_filtered.csv", row.names=1)
# read annotation file
annotation <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
##  Assign proper column names
num_samples <- ncol(Raw_counts)  # Get number of samples
colnames(Raw_counts) <- paste0("Sample_", seq_len(num_samples))  #  Raw_counts has no column names, this ensures they get unique names.


#  Add GeneID as a column for merging
Raw_counts$GeneID <- rownames(Raw_counts)

#  Merge with annotation to get Gene Symbols
merged_data <- merge(Raw_counts, annotation[, c("GeneID", "Symbol")], by = "GeneID", all.x = TRUE)

# Remove rows where Symbol is missing
merged_data_clean <- merged_data[!is.na(merged_data$Symbol), ]

#  Handle duplicate gene symbols by summing expression values
final_counts_agg <- merged_data_clean %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames(var = "Symbol")

# Convert final_counts_agg to a numeric matrix
final_counts_matrix <- as.matrix(final_counts_agg)
mode(final_counts_matrix) <- "numeric"

# Save final processed data
write.csv(final_counts_agg, "Processed_GSE224377.csv", row.names = TRUE, quote = FALSE)



#  Visualize distribution (Histogram of Log2-transformed counts)
hist(log2(as.numeric(unlist(final_counts_agg)) + 1), 
     col = "orange", 
     main = "Histogram of Processed Counts (log2)")

processed_data <- read.csv("Processed_GSE224377.csv", row.names=1)


#explore if is there any missing expression value (empty cell)
sum(is.na(processed_data))
is.na(processed_data)
is.null(processed_data)

#  Create Mapping from Sample_X to sample IDs


# Get the GSM IDs from the pheno file
gsm_ids <- rownames(pheno)

# Map Sample_1, Sample_2, ... to GSM7021440, GSM7021441, ...
# Ensure lengths match
if (length(gsm_ids) == ncol(processed_data)) {
  colnames(processed_data) <- gsm_ids
} else {
  stop("Number of samples in processed_data and pheno file do not match!")
}

# Verify the Alignment
# -----------------------------------------

# Check if sample names now match between processed_data and pheno for deseq2 analyssi
all(colnames(processed_data) %in% rownames(pheno))  # already return TRUE

################# DO the differential EXP analysis using DeSeq2##########
library(DESeq2)


#  Define Conditions

cond1 <- "lesion"  # Disease
cond2 <- "NAWM"    # Healthy

#  Create DESeq2 Dataset


# Ensure 'condition' and 'Patient.no' are factors
pheno$condition <- factor(pheno$condition)
pheno$Patient.no <- factor(pheno$Patient.no)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = processed_data,
                              colData = pheno,
                              design = ~ condition + Patient.no)

#  Run DESeq2 Workflow

dds <- DESeq(dds)


#  Extract Results (Specify Contrast)
# Contrast between lesion vs NAWM
res <- results(dds, contrast = c("condition", cond1, cond2))

#  Explore and Save Results
# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 1
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1,]
names_deges=rownames(deseq.deg)
#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)



# -----------------------------------------
#  Visualization 


#drow DEGs volcano plot
# Plot basic volcano
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="lesion vs NAWM DEGs", xlab="log2FoldChange", ylab="-log10(padj)"))

# Highlight upregulated genes (blue)
with(subset(res, padj < 0.05 & log2FoldChange > 1),
     points(log2FoldChange, -log10(padj), pch=20, col="blue"))

# Highlight downregulated genes (red)
with(subset(res, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(padj), pch=20, col="red"))

# legend 
legend("topright", legend = c("Upregulated", "Downregulated"), 
       col = c("blue", "red"), pch = 19, cex = 0.8, bty = "n")


#  Add Gene Names for Significant DEGs

# Define thresholds for labeling
sig_up <- subset(res, padj < 0.01 & log2FoldChange > 1)
sig_down <- subset(res, padj < 0.01 & log2FoldChange < -1)

# Add labels for upregulated genes
with(sig_up, text(log2FoldChange, -log10(padj), labels = rownames(sig_up), pos = 3, cex = 0.6, col = "blue"))

# Add labels for downregulated genes
with(sig_down, text(log2FoldChange, -log10(padj), labels = rownames(sig_down), pos = 3, cex = 0.6, col = "red"))
###############################################
##AHEATMAP
# FIRST  Normalize the data to make heatmap
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized = TRUE))

#  Extract DEGs expression values
exp.degs <- as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])

#  Prepare sample annotations
condition_ann <- data.frame(Condition = pheno$condition)
rownames(condition_ann) <- colnames(exp.degs)

#  Plot aHeatmap
aheatmap(log2(exp.degs + 1), 
         annCol = condition_ann, 
         main = "mRNA lesion vs NAWM", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         fontsize = 10)










