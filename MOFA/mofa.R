

# Load BiocManager
library(BiocManager)

# Install MOFA2 package
BiocManager::install("MOFA2", force = TRUE)

# Install devtools package
install.packages("devtools")


install.packages("ggplot2")
install.packages("randomForest")
install.packages("survminer")
# Load required libraries
library(devtools)
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reticulate)
library(randomForest)
library(utils)
library(survival)
library(survminer)
library(readxl)




# Load Data

epigenomic_data <- read.csv("DMPs_matrix.csv", row.names = 1)
transcriptomic_data <- read.csv("NormalizedCounts_Annotated.csv", row.names = 1)
colnames(transcriptomic_data)
# Load metadata for both datasets
Metad <-read_excel("metadata1.xlsx")
Metad$Sample_id_mrna<-gsub("-",".",Metad$Sample_id_mrna)
Metad$Sample_id_cpg <- gsub("-",".",Metad$Sample_id_cpg)## Replacing hyphens (-) with periods (.) makes the daata abvious and subsequent analysis.

common_samples <- intersect(colnames(transcriptomic_data), colnames(epigenomic_data))
mRNA_counts <- transcriptomic_data[, common_samples]
Methylation_counts <- epigenomic_data[, common_samples]
# Store in a named list for MOFA
all_omics_list <- list(
  mRNA = as.matrix(transcriptomic_data),
  Methylation = as.matrix(epigenomic_data)
)

# Ensure metadata sample names match data sample names
Metad <- Metad %>%
  filter(sample %in% common_samples)


colnames(transcriptomic_data)
Metad$ID
metaord= Metad[match(colnames(transcriptomic_data),Metad$Sample_id_mrna),]
colnames(transcriptomic_data) <- metaord$ID
transcriptomic_data <- transcriptomic_data[, !is.na(colnames(transcriptomic_data))] ##to remove the last na
metaord2= Metad[match(colnames(epigenomic_data),Metad$Sample_id_cpg),]
colnames(epigenomic_data) <- metaord2$ID

colnames(epigenomic_data)

# creating all omics list
all_omics_list <- list(data.matrix(epigenomic_data),data.matrix(transcriptomic_data))
names(all_omics_list)<-c("Methylation", "mRNA")


# Creating Mofa object ----
MOFAobject<- create_mofa(all_omics_list)

colnames(Metad)[7] <- "sample"
samples_metadata(MOFAobject) <- Metad
plot_data_overview(MOFAobject)

#### Fit MOFA model ----
## Defining Data Options 

data_opt <- get_default_data_options(MOFAobject)
data_opt
## Defining Model Options 

model_opt <- get_default_model_options(MOFAobject)
model_opt

## Defining Training Options 

train_opt <- get_default_training_options(MOFAobject)
train_opt$verbose <- T   ##to get output during analysis
train_opt
## Preparing MOFA
MOFAobject_prepared <- prepare_mofa(object = MOFAobject,
                                    data_options = data_opt,
                                    model_options = model_opt,
                                    training_options = train_opt)


## Running MOFA
MOFAobject_run <- run_mofa(MOFAobject_prepared,use_basilisk = T)
##succsfully converged at iteration 11 according to elbo



# Save MOFA model
saveRDS(MOFAobject, "MOFAobject.rds")

#  Visualize Results --------

# Load trained MOFA model
MOFAobject <- readRDS("MOFAobject.rds")

# Extract factors
df_factors <- get_factors(MOFAobject, factors = "all", as.data.frame = TRUE)

# Merge with metadata
df_plot <- left_join(df_factors, Metad, by = c("Sample" = "sample"))

# Plot first two factors
ggplot(df_plot, aes(x = Factor1, y = Factor2, color = Sample_source_name_ch1_cpg)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("MOFA Integration of mRNA and Methylation Data")

# Identifying the top factors----
plot_factor_cor(MOFAobject_run) ##check model quality

correlate_factors_with_covariates(MOFAobject_run, 
                                  covariates = c("Sample_id_cpg",
                                                 "Sample_id_mrna",
                                                 "Sample_title_cpg",
                                                 "Sample_source_name_ch1_cpg",
                                                 "age",
                                                 "gender",
                                                 "sample"
                                  ), 
                                  plot="log_pval") #means the plot will display log-transformed p-values which higher p value is more signficant
# Sample_source name
plot_factor(
  MOFAobject_run,
  factors = 1,
  color_by   = "Sample_source_name_ch1_cpg",
  dot_size = 3,
  dodge = T,  #Ensures the dots are slightly separated also avoiding overlape
  add_boxplot   = T
)

plot_factor(
  MOFAobject_run,
  factors = 2,
  color_by   = "Sample_source_name_ch1_cpg",
  dot_size = 3,
  dodge = T,  #Ensures the dots are slightly separated also avoiding overlape
  add_boxplot   = T
)

#gender
plot_factor(
  MOFAobject_run,
  factors = 3,
  color_by   = "gender",
  dot_size = 3,
  dodge = T,  #Ensures the dots are slightly separated also avoiding overlape
  add_boxplot   = T
)

#sample title cpg
plot_factor(
  MOFAobject_run,
  factors = 3,
  color_by   = "Sample_title_cpg",
  dot_size = 3,
  dodge = T,  #Ensures the dots are slightly separated also avoiding overlape
  add_boxplot   = T
)
#sample

plot_factor(
  MOFAobject_run,
  factors = 3,
  color_by   = "sample",
  dot_size = 3,
  dodge = T,  #Ensures the dots are slightly separated also avoiding overlape
  add_boxplot   = T
)

#identifying the different omics weights----

factor_var <- get_variance_explained(MOFAobject_run) 
factor_var$r2_per_factor
plot_variance_explained(MOFAobject_run)
plot_variance_explained(MOFAobject_run,plot_total = T) 


##now w make heatmap
##Understanding relationships between factors and samples.,Identifying clusters or patterns in the dataset.
plot_weights_heatmap(MOFAobject_run,
                     view = "Methylation",
                     show_colnames=F)



plot_variance_explained(MOFAobject_run,factors = 1) #cpg
plot_variance_explained(MOFAobject_run,factors = 2)  #cpg
plot_variance_explained(MOFAobject_run,factors = 3)  #cpg

# identifying the different features weights----
# feature weights for factor 1
plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 1, 
  nfeatures = 10
)

## ineed now to identfy each factor what is top features??
plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 1, 
  nfeatures = 10
)
plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 1,  
                  features = 10,
                  sign = "all",
                  color_by = "Sample_source_name_ch1_cpg"
) + labs(y="mRNA")




# identifying the different features weights----
# feature weights for factor 1
plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 2, 
  nfeatures = 10
)

## ineed now to identfy each factor what is top features??
plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 2, 
  nfeatures = 10
)
plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 2,  
                  features = 10,
                  sign = "all",
                  color_by = "Sample_source_name_ch1_cpg"
) + labs(y="mRNA")


plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 3, 
  nfeatures = 10
)

## ineed now to identfy each factor what is top features??
plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 3, 
  nfeatures = 10
)
plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 3,  
                  features = 10,
                  sign = "all",
                  color_by = "Sample_source_name_ch1_cpg"
) + labs(y="mRNA")





#now make heatmap to  validate how Factor 1 (and the topfeatures) aligns with Sample_source_name_ch1_cpg
# HEATMAPS----

plot_data_heatmap(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 1, 
  features = 10, 
  show_rownnames = T,
  show_colnames=F,
  denoise = TRUE, 
  annotation_samples = "Sample_source_name_ch1_cpg",
  show_colnnames= F
)

#########
plot_top_weights(
  MOFAobject_run, 
  view = "mRNA", 
  factor = 1, 
  nfeatures = 10
)
plot_data_scatter(MOFAobject_run, 
                  view = "mRNA",
                  factor = 1,  
                  features = 10,
                  sign = "all",
                  color_by = "Sample_source_name_ch1_cpg"
) + labs(y="mRNA")
###########################
