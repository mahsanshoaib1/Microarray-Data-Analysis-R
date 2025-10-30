# end to end microarray gene expression analysis project with complete details and notes
# Grade 0

#---------------------------install all the packages required-----------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")        # Helps install bioinformatics packages
# install.packages("affy")       # to read raw microarray data (CEL files)
# install.packages("GEOquery")   # To fetch public datasets from GEO database
# install.packages("limma")      # For differential expression analysis
# install.packages("Biobase")
# install.packages("oligo")
# install.packages("pheatmap")
# install.packages("stringr")
# install.packages("pd.hta.2.0")
# install.packages("ggplot2")
# install.packages("annotate")     # To map gene names. Helps in converting probe IDs to gene symbols
## if the installation give error use the BiocManager as given below
# BiocManager::install("affy")
# if (!requireNamespace("hta20transcriptcluster.db", quietly = TRUE))
#   BiocManager::install("hta20transcriptcluster.db")


#---------------------------Activate all the packages----------------------------

library(ggplot2) # to draw plots
library(limma) # for microarray analysis
library(Biobase)
library(BiocManager)
library(affy)
library(oligoClasses)
library(oligo)
library(pheatmap)
library(stringr)
library(pd.hta.2.0)
library(GEOquery)
library(annotate)
library(data.table)
library(tidyverse)


#---------------------------Setting Directory and file paths------------------

path <- setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_2")
setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_2")


#---------------------------Playing with metadata/sdrf dataset------------------

sdrf <- read.delim("E-GEOD-76250.sdrf.txt", header = TRUE, sep = "\t")

cols <- colnames(sdrf) # to get the names of columns
nrow(sdrf) # to get the number of rows

# extracting rows in which tumor grade is equal to 0 and tissue type if normal and saving it in a variable.
Grade_0 <- filter(sdrf, FactorValue..tumor.grade.== 0 | Characteristics..tissue.type.== "normal breast tissue")
unique(Grade_0$FactorValue..tumor.grade.) # checking unique values in the tumor grade columns
unique(Grade_0$Array.Data.File) 
count(Grade_0)
head(Grade_0)


#---------------------------Loading and reading the .CEL dataset------------------------#

# rownames(Grade_0) <- Grade_0$Array.Data.File # making Array.Data.File column as rownames


# Reading CEL FILES using Grade_0 data
Grade_0 <- AnnotatedDataFrame(Grade_0)
local_data <-oligo::read.celfiles(filenames = file.path(path,Grade_0$Array.Data.File), phenoData = Grade_0)
stopifnot(validObject(local_data))
colnames(pData(local_data))
rownames(pData(local_data))
# rownames(pData(local_data)) <- Grade_0$Array.Data.File # making Array.Data.File column as rownames
# rownames(pData(local_data))


# taking important columns of metadata
head(Biobase::pData(local_data))
Biobase::pData(local_data) <- Biobase::pData(local_data)[,
                                                         c("Source.Name", "Array.Data.File", "FactorValue..tumor.grade.", "Characteristics..tumor.grade.", "Characteristics..tissue.type.")]

colnames(pData(local_data))


#---------------------------Checking the Quality of RAW Data-----------------------#

## boxplot  
# ?boxplot
oligo::boxplot(local_data, target = "core", main="Grade 0 TNBC Raw Data")

## checking null values
sum(is.na(exprs(local_data)))

## summary of data
raw_data_summary <- summary(exprs(local_data)) # giving the error due to size of data
# view(data_summary)

setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/grade_0")

write.csv(raw_data_summary, "data_summary_Grade_0_tnbc_raw.csv")


## histogram: Helps identify batch effects or abnormal distributions.
# oligo::hist(exprs(local_data))
hist(as.vector(exprs(local_data)), breaks = 100, main = "Intensity Distribution of Grade 0 Raw Data") 

## density plot
affy::hist(local_data, target= "core", col = 1:length(sampleNames(local_data)), lty = 1, main = "Density Plot of Grade 0 TNBC Raw Data")

## log transformation of the data
# exp_raw <- log2(exprs(local_data)) # rma autometically perform log2 transformation
# view(exp_raw)

## MA plot: mean vs log fold change, global intensity differences between arrays.
# MAplot(local_data)

## RNA degradation plot
# deg <- AffyRNAdeg(local_data)
# plotAffyRNAdeg(deg)  #High degradation suggests poor RNA sample quality.

## heatmap
pheatmap(cor(exprs(local_data)))
## Hierarchical clustering
hc <- hclust(dist(t(exprs(local_data))))
plot(hc, main = "Hierarchical Clustering Dendrogram of Grade 0 TNBC Raw Data")


#---------------------------Normalization---------------------------
  
## Normalize using RMA
norm_data <- rma(local_data)
expr_matrix <- exprs(norm_data)

#---------------------------Checking the Quality of Normalized Data---------------------------

# box plot
oligo::boxplot(expr_matrix, target = "core", main="Grade 0 TNBC Normalized Data Quality Check")

## checking null values
sum(is.na(expr_matrix))

# summary
normalized_data_summary <- summary(expr_matrix) # giving the error due to size of data
write.csv(normalized_data_summary, "normalized_data_summary_Grade_0_tnbc_raw.csv")

# histogram: Helps identify batch effects or abnormal distributions.
# oligo::hist(exprs(norm_data))
hist(as.vector(exprs(norm_data)), breaks = 100, main = "Intensity Distribution of Grade 0 TNBC Normalized Data") # giving the error due to size of data

#density plot
affy::hist(norm_data, col = 1:length(sampleNames(norm_data)), lty = 1, main = "Density Plot of Grade 0 TNBC Normalized Data")

# Principle Component Analysis
pca <- prcomp(t(expr_matrix), scale. = TRUE)

percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2), 1) # Calculate Percentage Variance for Each Principal Component
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     phenotype = pData(norm_data)$FactorValue..tumor.grade.) # this code is creating the dataframe of pca to plot it in the next step. 
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = phenotype, color = phenotype)) +
  ggtitle("PCA plot of the Grade 0 TNBC Normalized epxression data")+
  xlab(paste0("PC1, varExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, varExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4, 15, 16, 17)) +
  scale_color_manual(values = c("red", "blue", "green", "purple"))
# the PCA part isnt running because of following error: Error: cannot allocate vector of size 6.5 Gb

# MA plot: mean vs log fold change, global intensity differences between arrays.
# MAplot(norm_data)

# RNA degradation plot
# deg <- AffyRNAdeg(norm_data)
# plotAffyRNAdeg(deg)  #High degradation suggests poor RNA sample quality.

# heatmap
pheatmap(cor(expr_matrix))

# Hierarchical clustering
hc <- hclust(dist(t(expr_matrix)))
plot(hc, main = "Hierarchical Clustering Dendrogram of Grade 0 TNBC Normalized Data")

# chatgpt used for this

# Remaining Steps in DGEA
# Annotation of Probes – Convert probe IDs to gene symbols.
# Filtering of Low-Quality Probes – Remove non-informative or low-expression probes.
# Setting Up the Experimental Design – Define groups (e.g., tumor vs. normal).
# Differential Expression Analysis – Use limma to find significantly differentially expressed genes (DEGs).
# Adjusting for Multiple Comparisons – Use False Discovery Rate (FDR) correction.
# Visualization of DEGs – Volcano plot, heatmaps.
# Functional Enrichment Analysis – Find pathways affected by DEGs.

#--------------------------- Step 1: Annotation (Convert Probe IDs to Gene Symbols) ---------------------------#

# Load annotation package for HTA 2.0

# Probe Annotation Practice. 
library(hta20transcriptcluster.db)
keytypes(hta20transcriptcluster.db)
columns(hta20transcriptcluster.db)
ls("package:hta20transcriptcluster.db")

# ?hta20transcriptcluster.db

#---------------------------------------#
# Final annotation code

# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data) # Get probe IDs
head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
gene_entrezid <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]

# Create an annotated expression matrix
expr_matrix <- exprs(norm_data) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

# view(expr_matrix)

# Check mapping for DDX11L1
# mapIds(hta20transcriptcluster.db, keys = "DDX11L1", column = c("ENTREZID"), keytype = "SYMBOL")
# mapIds(hta20transcriptcluster.db, keys = "LINC01001", column = c("ENTREZID"), keytype = "SYMBOL")
# mapIds(hta20transcriptcluster.db, keys = "OR4F5", column = c("ENTREZID"), keytype = "SYMBOL")

# Convert to a data frame (required to manipulate column names easily)
type(expr_matrix)
nrow(expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix %>% distinct(SYMBOL, .keep_all = TRUE)
view(head(expr_matrix))
type(expr_matrix)
nrow(expr_matrix)

# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
view(head(expr_matrix))

# write.csv(expr_matrix, "final_Grade_0_TNBC_Normalized_Data.csv")


#--------------------------- Step 2: Filter Low-Quality Probes ---------------------------#

# Remove low-expressed genes (genes with very low variance)

filtered_matrix <- expr_matrix
# view(head(filtered_matrix))
nrow(filtered_matrix)

# Store the gene symbols (row names) before conversion
gene_symbols <- rownames(filtered_matrix)

# Convert all columns to numeric, preserving row names
filtered_matrix <- as.data.frame(lapply(filtered_matrix, as.numeric))
rownames(filtered_matrix) <- gene_symbols

# Verify the result
# view(head(filtered_matrix))

#--------------------------- Step 3: Define Experimental Design ---------------------------#
dim(filtered_matrix)  # Check dimensions (rows = genes, columns = samples)
dim(design)           # Check dimensions (rows = experimental conditions, columns = factors)

# Extract phenotype data from `sdrf`
group <- factor(Grade_0$FactorValue..tissue.type.) # Normal vs Tumor
design <- model.matrix(~ 0 + group)
levels(group)
colnames(design) <- levels(group)
rownames(design) <- colnames(filtered_matrix)

#--------------------------- Step 4: Differential Expression Analysis (Using Limma) ---------------------------#

library(limma)

# Create a contrast matrix (comparing Tumor vs Normal)
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(grade0TNBC_vs_Normal = `triple.negative.breast.cancer..TNBC.` - `normal.breast.tissue`, levels = design)

class(filtered_matrix)
# Fit linear model
fit <- lmFit(filtered_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get top differentially expressed genes
deg_results <- topTable(fit2, coef = "grade0TNBC_vs_Normal", number = Inf, adjust.method = "fdr")

# Save the results
# write.csv(deg_results, "Grade_0_TNBC_DEG_results.csv", row.names = TRUE)

dim(expr_matrix)  # Should show (genes, samples)
head(expr_matrix)

#most significant DEG
summary(deg_results$adj.P.Val)

DEG_norm_DS <- subset(deg_results, adj.P.Val < 0.05)

nrow(DEG_norm_DS)
view(head(DEG_norm_DS))
write.csv(DEG_norm_DS, "updated_Grade_0_TNBC_DEG_results.csv")

#--------------------------- Step 5: Visualization ---------------------------#

#Now time to generate a volcano plot
str(deg_results$P.Value)
colnames(deg_results)
sum(is.na(deg_results$P.Value))

deg_results$P.Value <- as.numeric(deg_results$P.Value)


## 1. Volcano Plot
library(ggplot2)

deg_results$Significance <- ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, "Significant", "Not Significant")

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("red", "gray")) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes of Grade 0 TNBC") +
  xlab("Log Fold Change") +
  ylab("-log10 Adjusted P-value")

## 2. Heatmap of Top 50 DEGs
library(pheatmap)

top50_genes <- rownames(deg_results)[1:50]
# pheatmap(filtered_matrix[top50_genes, ], scale = "row", annotation_col = data.frame(Group = group)) # giving error



###########################################################################

#-----------------------------------------------------------------------#
# Heat map

# Updated Heatmap with Euclidean Distance and Better Colors
library(pheatmap)
library(stringr)

colnames(pData(norm_data))

disease_names <- ifelse(str_detect(pData(norm_data)$Characteristics..tissue.type.,
                                   "normal.breast.tissue"), "normal.breast.tissue", "triple.negative.breast.cancer..TNBC.")

annotation_heatmap <- data.frame(disease= disease_names)
annotation_heatmap
row.names(annotation_heatmap) <-row.names(pData(norm_data))

row.names(pData(norm_data))
row.names(annotation_heatmap)
# Recalculate distances using Euclidean
dists <- as.matrix(dist(t(expr_matrix), method = "euclidean"))
dists
rownames(dists) <- row.names(pData(norm_data))
rownames(dists)

# Adjust color palette for better contrast
hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(225)

colnames(dists) <- NULL
diag(dists) <- NA


# Define annotation colors for stages
#Updated Heatmap with Euclidean Distance and Better Colors
library(pheatmap)
library(stringr)

pData(norm_data)
disease_names <- case_when(
  str_detect(pData(norm_data)$Characteristics..tissue.type., "normal.breast.tissue") ~ "normal",
  str_detect(pData(norm_data)$Characteristics..tissue.type., "triple.negative.breast.cancer..TNBC.") ~ "TNBC"
)


annotation_heatmap <- data.frame(disease= disease_names)
annotation_heatmap
row.names(annotation_heatmap) <-row.names(pData(norm_data))

row.names(pData(norm_data))
row.names(annotation_heatmap)
# Recalculate distances using Euclidean
dists <- as.matrix(dist(t(expr_matrix), method = "euclidean"))
dists
rownames(dists) <- row.names(pData(norm_data))
rownames(dists)

# Adjust color palette for better contrast
hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(225)


colnames(dists) <- NULL
diag(dists) <- NA


# Define annotation colors for stages
ann_colors <- list(
  disease_name = c(
    normal = "green", TNBC = "orange"
  )
)

# Generate improved heatmap
pheatmap(dists, col = hmcol,
         annotation_row = annotation_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 30,
         cluster_rows = TRUE,  # Ensure proper clustering
         cluster_cols = TRUE,
         legend_breaks = c(min(dists, na.rm = TRUE), max(dists, na.rm = TRUE)),
         legend_labels = c("small distance", "large distance"),
         main = "Improved Clustering Heatmap of Grade 0 TNBC Differentially expressed genes")




##############################################################################
##############################################################################
##############################################################################
