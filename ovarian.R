# end to end microarray gene expression analysis project with complete details and notes
# Ovarian Cancer


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
# if (!requireNamespace("hgu133a.db", quietly = TRUE))
#   BiocManager::install("hgu133a.db")


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
library(pd.hg.u133a)
library(GEOquery)
library(annotate)
library(data.table)
library(tidyverse)


#---------------------------Setting Directory and file paths------------------

path <- setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_6")
setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_6")


#---------------------------Playing with metadata/ovarian dataset------------------

ovarian <- read.delim("E-GEOD-26712.sdrf.txt", header = TRUE, sep = "\t")

cols <- colnames(ovarian) # to get the names of columns
nrow(ovarian) # to get the number of rows

#---------------------------Loading and reading the .CEL dataset------------------------#

# reading cel files
local_data <-oligo::read.celfiles(list.celfiles()) 
stopifnot(validObject(local_data))

# Extract filenames from local_data
local_files <- basename(colnames(local_data))  # Extract just file names, ignoring paths
local_files

# Filter ovarian data to keep only rows where "Array Data File" exists in local_files
ovarian_filtered <- ovarian[ovarian$`Array.Data.File` %in% local_files, ]
head(ovarian_filtered)
nrow(ovarian_filtered)
write.table(ovarian_filtered, "filtered.sdrf.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ovarian <- read.delim("filtered.sdrf.txt", header = TRUE, sep = "\t")

# Ensure column names match correctly
colnames(ovarian) <- make.names(colnames(ovarian))  # Convert column names to valid R names

# Ensure that sample names in `ovarian` match `local_data`
rownames(ovarian)
rownames(ovarian) <- basename(ovarian$`Array.Data.File`)  # Set row names from file names
rownames(ovarian)

# Assign phenotype data to the existing local_data object
pData(local_data) <- ovarian

# taking important columns of metadata
head(Biobase::pData(local_data))
Biobase::pData(local_data) <- Biobase::pData(local_data)[,
                                                         c("Source.Name", "Array.Data.File", "Characteristics.tissue.")]

colnames(pData(local_data))

local_data
#---------------------------Checking the Quality of RAW Data-----------------------#

## boxplot
oligo::boxplot(local_data, target = "core", main="ovarian Raw Data")

setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/ovarian")

## checking null values
sum(is.na(exprs(local_data)))

## summary of data
#raw_data_summary <- summary(exprs(local_data)) # giving the error due to size of data
#write.csv(raw_data_summary, "data_summary_ovarian_raw.csv")


## histogram: Helps identify batch effects or abnormal distributions.
# oligo::hist(exprs(local_data))
hist(as.vector(exprs(local_data)), breaks = 100, main = "Intensity Distribution of ovarian Raw Data") 

## density plot
affy::hist(local_data, target= "core", col = 1:length(sampleNames(local_data)), lty = 1, main = "Density Plot of ovarian  Raw Data")

## log transformation of the data
# exp_raw <- log2(exprs(local_data)) # rma autometically perform log2 transformation
# view(exp_raw)

## MA plot: mean vs log fold change, global intensity differences between arrays.
# MAplot(local_data)

## RNA degradation plot
# deg <- AffyRNAdeg(local_data)
# plotAffyRNAdeg(deg)  #High degradation suggests poor RNA sample quality.

## heatmap
# pheatmap(cor(exprs(local_data)))
# ## Hierarchical clustering
# hc <- hclust(dist(t(exprs(local_data))))
# plot(hc, main = "Hierarchical Clustering Dendrogram of ovarian  Raw Data")


#---------------------------Normalization---------------------------

## Normalize using RMA
norm_data <- rma(local_data)
view(head(norm_data))
expr_matrix <- exprs(norm_data)
view(head(expr_matrix))
#---------------------------Checking the Quality of Normalized Data---------------------------

# box plot
oligo::boxplot(expr_matrix, target = "core", main="ovarian  Normalized Data Quality Check")

## checking null values
# sum(is.na(expr_matrix))

# summary
#normalized_data_summary <- summary(expr_matrix) # giving the error due to size of data
#write.csv(normalized_data_summary, "normalized_data_summary_ovarian_raw.csv")

# histogram: Helps identify batch effects or abnormal distributions.
# oligo::hist(exprs(norm_data))
hist(as.vector(exprs(norm_data)), breaks = 100, main = "Intensity Distribution of ovarian  Normalized Data") # giving the error due to size of data

#density plot
affy::hist(norm_data, col = 1:length(sampleNames(norm_data)), lty = 1, main = "Density Plot of ovarian  Normalized Data")

# Principle Component Analysis
# pca <- prcomp(t(expr_matrix), scale. = TRUE)
# 
# percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2), 1) # Calculate Percentage Variance for Each Principal Component
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])
# dataGG <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
#                      phenotype = pData(norm_data)$Characteristics.tissue.) # this code is creating the dataframe of pca to plot it in the next step. 
# ggplot(dataGG, aes(PC1, PC2)) +
#   geom_point(aes(shape = phenotype, color = phenotype)) +
#   ggtitle("PCA plot of the ovarian  Normalized epxression data")+
#   xlab(paste0("PC1, varExp: ", percentVar[1], "%")) +
#   ylab(paste0("PC2, varExp: ", percentVar[2], "%")) +
#   theme(plot.title = element_text(hjust = 0.5))+
#   coord_fixed(ratio = sd_ratio) +
#   scale_shape_manual(values = c(4, 15, 16, 17)) +
#   scale_color_manual(values = c("red", "blue", "green", "purple"))
# the PCA part isnt running because of following error: Error: cannot allocate vector of size 6.5 Gb

# MA plot: mean vs log fold change, global intensity differences between arrays.
# MAplot(norm_data)

# RNA degradation plot
# deg <- AffyRNAdeg(norm_data)
# plotAffyRNAdeg(deg)  #High degradation suggests poor RNA sample quality.

# heatmap
#pheatmap(cor(expr_matrix))

# Hierarchical clustering
# hc <- hclust(dist(t(expr_matrix)))
# plot(hc, main = "Hierarchical Clustering Dendrogram of ovarian  Normalized Data")

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


# Probe Annotation Practice. 
library(hgu133a.db)
keytypes(hgu133a.db)
columns(hgu133a.db)
ls("package:hgu133a.db")

?hgu133a.db

#---------------------------------------#
library(hgu133a.db)
keytypes(hgu133a.db)
columns(hgu133a.db)
ls("package:hgu133a.db")


#---------------------------------------#
# Final annotation code

# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data) # Get probe IDs
head(keys(hgu133a.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hgu133a.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
gene_entrezid <- mapIds(hgu133a.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]
gene_symbols
# Create an annotated expression matrix
expr_matrix <- exprs(norm_data) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

view(head(expr_matrix))

# Check mapping for DDX11L1
mapIds(hgu133a.db, keys = "DDX11L1", column = c("ENTREZID"), keytype = "SYMBOL")
mapIds(hgu133a.db, keys = "LINC01001", column = c("ENTREZID"), keytype = "SYMBOL")
mapIds(hgu133a.db, keys = "OR4F5", column = c("ENTREZID"), keytype = "SYMBOL")

# Convert to a data frame (required to manipulate column names easily)
type(expr_matrix)
nrow(expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
view(head(expr_matrix))
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix[!duplicated(expr_matrix$SYMBOL), ]
view(head(expr_matrix))
type(expr_matrix)
nrow(expr_matrix)


# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
view(head(expr_matrix))

# write.csv(expr_matrix, "Grade_2_TNBC_Normalized_Data.csv")


#--------------------------- Step 2: Filter Low-Quality Probes ---------------------------#

# Remove low-expressed genes (genes with very low variance)

filtered_matrix <- expr_matrix
view(head(expr_matrix))
nrow(filtered_matrix)
view(head(expr_matrix))
# Store the gene symbols (row names) before conversion
gene_symbols <- rownames(filtered_matrix)
gene_symbols
# Convert all columns to numeric, preserving row names
filtered_matrix <- as.data.frame(lapply(filtered_matrix, as.numeric))
rownames(filtered_matrix) <- gene_symbols

# filtered_matrix <- expr_matrix
# # view(head(filtered_matrix))
# filtered_matrix <- as.data.frame(sapply(filtered_matrix, as.numeric))
# keep_genes <- rowMeans(filtered_matrix) > 5  # Adjust threshold as needed
# filtered_matrix <- filtered_matrix[keep_genes, ]
# filtered_matrix$ENTREZID <- NULL

#--------------------------- Step 3: Define Experimental Design ---------------------------#

filtered_matrix <- subset(filtered_matrix, select = -ENTREZID)
view(head(filtered_matrix))
# Extract phenotype data from `ovarian`
group <- factor(ovarian$Characteristics.tissue.) # Normal vs Tumor
design <- model.matrix(~ 0 + group)
levels(group)

dim(design)
dim(filtered_matrix)
colnames(design) <- levels(group)
rownames(design) <- colnames(filtered_matrix)
view(head(filtered_matrix))
#--------------------------- Step 4: Differential Expression Analysis (Using Limma) ---------------------------#

library(limma)

# Create a contrast matrix (comparing Tumor vs Normal)
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(ovarian_vs_Normal = `Late.stage.high.grade.ovarian.cancer` - `Normal.ovarian.surface.epithelium`, levels = design)

# Fit linear model
fit <- lmFit(filtered_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get top differentially expressed genes
deg_results <- topTable(fit2, coef = "ovarian_vs_Normal", number = Inf, adjust.method = "fdr")

# Save the results
# write.csv(deg_results, "ovarian_DEG_results.csv", row.names = TRUE)

dim(expr_matrix)  # Should show (genes, samples)
head(expr_matrix)

#most significant DEG
summary(deg_results$adj.P.Val)

DEG_norm_DS <- subset(deg_results, adj.P.Val < 0.05)

nrow(DEG_norm_DS)
view(head(DEG_norm_DS))
write.csv(DEG_norm_DS, "updated_ovarian_DEG_results.csv")

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
  ggtitle("Volcano Plot of Differentially Expressed Genes of ovarian ") +
  xlab("Log Fold Change") +
  ylab("-log10 Adjusted P-value")

## 2. Heatmap of Top 50 DEGs
library(pheatmap)

top50_genes <- rownames(deg_results)[1:50]
# pheatmap(filtered_matrix[top50_genes, ], scale = "row", annotation_col = data.frame(Group = group)) # giving error



###########################################################################

#-----------------------------------------------------------------------#
# Heat map

# Define annotation colors for stages
#Updated Heatmap with Euclidean Distance and Better Colors
library(pheatmap)
library(stringr)

pData(norm_data)
disease_names <- case_when(
  str_detect(pData(norm_data)$Characteristics.tissue., "Normal.ovarian.surface.epithelium") ~ "normal",
  str_detect(pData(norm_data)$Characteristics.tissue., "Late.stage.high.grade.ovarian.cancer") ~ ""
)


annotation_heatmap <- data.frame(disease= disease_names)
annotation_heatmap
row.names(annotation_heatmap) <-row.names(pData(norm_data))

row.names(pData(norm_data))
row.names(annotation_heatmap)

# Recalculate distances using Euclidean
dists <- as.matrix(dist(t(expr_matrix), method = "euclidean"))
rownames(dists) <- row.names(pData(norm_data))

# Adjust color palette for better contrast
hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(225)


colnames(dists) <- NULL
diag(dists) <- NA


# Define annotation colors for stages
ann_colors <- list(
  disease_name = c(
    normal = "green", moderate = "orange"
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
         main = "Improved Clustering Heatmap of ovarian  Differentially expressed genes")




##############################################################################
##############################################################################
##############################################################################
