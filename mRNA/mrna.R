
library(BiocManager)

library(DESeq2)

library(NMF)

#loading the gene expression data 


md <- as.matrix(read.xlsx("GSE224377_raw_counts_GRCh38.p13_NCBI.xlsx", sheet = 1, rowNames=TRUE)
duplicated(md$GeneID)         
#loading the phenotype data 

pheno <- read.xlsx("pheno.xlsx", sheet = 1,rowNames = TRUE)
table(pheno$condition)

#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(data)
dim(pheno)
#explore the data distribution using the histogram plot
hist(data, col = "orange", main="Histogram")

#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(data+1), col = "orange", main="Histogram")

boxplot(log2(data[1:5,]+1))

# QQ plot for the normality
qqnorm(data[3,])
qqline(data[3,])

shapiro.test(data[100,])

#explore if is there any missing expression value (empty cell)
sum(is.na(data))
is.na(data)
is.null(data)
sum(is.nan(data))

#It is absolutely critical that the columns of the "data" and the rows of 

###### DO the differential EXP analysis using DeSeq2

#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="lesion" #disease
cond2="NAWM"    #healthy

#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = data , colData = pheno, design = ~ condition +Patient.no)

#run the deseq2 worflow
dds.run = DESeq(dds)
#specifying teh contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("condition",cond1 ,cond2))

# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1,]
names_deges=rownames(deseq.deg)
#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)

#drow DEGs volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="lesion vs NAWM DEGs"))
with(subset(res, padj<.05 & (log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-3.7,y=2.5,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

degs_normalizedcount=normalized_counts[names_deges,]
write.csv(degs_normalizedcount, file = "DEGs_NormalizedCounts.csv", row.names = TRUE, quote = FALSE)

####drow heatmap####
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
# Save to a CSV file
write.csv(normalized_counts, "normalized_counts.csv", row.names=TRUE)
#extract counts values of DEGs only for each stage
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
heatmap(log2(exp.degs+1))
#another advanced heatmap
aheatmap(log2(exp.degs+1), annCol =pheno$condition, main="mRNA lesion vs NAWM")


####reading expression annotation file 
library(readr)
annotatgenes=as.data.frame(read_tsv("Human.GRCh38.p13.annot.tsv"))
# Aggregate by mean 

annotatgenes <- annotatgenes[!duplicated(annotatgenes$Symbol), ]
rownames(normalized_counts) <- annotatgenes$Symbol
annotatgenes$Symbol <- make.names(annotatgenes$Symbol, unique = TRUE)
rownames(normalized_counts) <- annotatgenes$Symbol
nrow(annotatgenes)  # Number of rows in annotatgenes
nrow(normalized_counts)  # Number of rows in normalized_counts
# Find indices of extra rows in normalized_counts
extra_rows <- setdiff(rownames(normalized_counts), rownames(annotatgenes))
length(extra_rows)  # Should show how many rows are mismatched
# Keep only rows in normalized_counts that match annotatgenes
normalized_counts <- normalized_counts[1:nrow(annotatgenes), ]
# Add Symbol column for aggregation
normalized_counts$Symbol <- annotatgenes$Symbol

# Aggregate to handle duplicates
aggregated_counts <- aggregate(. ~ Symbol, data = normalized_counts, FUN = mean)

# Set row names and remove Symbol column
rownames(aggregated_counts) <- aggregated_counts$Symbol
aggregated_counts <- aggregated_counts[, -1]

# View the final result
head(aggregated_counts)


annotatgenes$Symbol[duplicated(annotatgenes$Symbol)]
rownames(normalized_counts) <- annotatgenes$Symbol

write.csv(normalized_counts, file = "NormalizedCounts_Annotated.csv", row.names = TRUE, quote = FALSE)




