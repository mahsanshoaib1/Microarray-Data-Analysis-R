# undestading microarray data
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



path <- setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_2")
# setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/data_2")


#---------------------------loading data------------------

sdrf <- read.delim("E-GEOD-76250.sdrf.txt", header = TRUE, sep = "\t")

# cols <- colnames(sdrf) # to get the names of columns
nrow(sdrf) # to get the number of rows

# extracting rows in which tumor grade is equal to 0 and tissue type if normal and saving it in a variable.
Grade_0 <- filter(sdrf, FactorValue..tumor.grade.== 0)
Grade_2 <- filter(sdrf, FactorValue..tumor.grade.== 2)
Grade_3 <- filter(sdrf, FactorValue..tumor.grade.== 3)
normal <- filter(sdrf, Characteristics..tissue.type.== "normal breast tissue")


#---------------------------Loading and reading the .CEL dataset---------------------------#

# rownames(Grade_0) <- Grade_0$Array.Data.File # making Array.Data.File column as rownames


# Reading CEL FILES using Grade_0 data
Grade_0 <- AnnotatedDataFrame(Grade_0)
Grade_2 <- AnnotatedDataFrame(Grade_2)
Grade_3 <- AnnotatedDataFrame(Grade_3)
Grade_normal <- AnnotatedDataFrame(normal)

data_0 <-oligo::read.celfiles(filenames = file.path(path,Grade_0$Array.Data.File), phenoData = Grade_0)
stopifnot(validObject(data_0))

data_2 <-oligo::read.celfiles(filenames = file.path(path,Grade_2$Array.Data.File), phenoData = Grade_2)
stopifnot(validObject(data_2))

data_3 <-oligo::read.celfiles(filenames = file.path(path,Grade_3$Array.Data.File), phenoData = Grade_3)
stopifnot(validObject(data_3))

data_normal <-oligo::read.celfiles(filenames = file.path(path,Grade_normal$Array.Data.File), phenoData = Grade_normal)
stopifnot(validObject(data_normal))


#---------------------------Probe Annotation of grade 0---------------------------#

# Probe Annotation Practice. 
library(hta20transcriptcluster.db)
keytypes(hta20transcriptcluster.db)
columns(hta20transcriptcluster.db)
ls("package:hta20transcriptcluster.db")


#---------------------------------------#
# Final annotation code of grade 0
norm_data_0 <- rma(data_0)
expr_matrix <- exprs(norm_data_0)


# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data_0)

head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
# gene_entrezid <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
# gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]

# Create an annotated expression matrix
expr_matrix <- exprs(norm_data_0) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

# Convert to a data frame (required to manipulate column names easily)
#type(expr_matrix)
#nrow(expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
#nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix %>% distinct(SYMBOL, .keep_all = TRUE)

type(expr_matrix)
nrow(expr_matrix)

# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
expr_matrix_0 <- expr_matrix
expr_matrix_0 <- t(expr_matrix_0)
view(head(expr_matrix_0))

#---------------------------Probe Annotation of grade 2---------------------------#

# Final annotation code of grade 2
norm_data_2 <- rma(data_2)
expr_matrix <- exprs(norm_data_2)


# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data_2)

head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
# gene_entrezid <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
# gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]

# Create an annotated expression matrix
expr_matrix <- exprs(norm_data_2) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

# Convert to a data frame (required to manipulate column names easily)
type(expr_matrix)
nrow(expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix %>% distinct(SYMBOL, .keep_all = TRUE)

# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
expr_matrix_2 <- expr_matrix
expr_matrix_2 <- t(expr_matrix_2)
#view(head(expr_matrix_2))

#---------------------------------------#
                        # Final annotation code of grade 3
norm_data_3 <- rma(data_3)
expr_matrix <- exprs(norm_data_3)


# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data_3)

head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
# gene_entrezid <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
# gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]

# Create an annotated expression matrix
expr_matrix <- exprs(norm_data_3) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

# Convert to a data frame (required to manipulate column names easily)
type(expr_matrix)
nrow(expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix %>% distinct(SYMBOL, .keep_all = TRUE)

# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
expr_matrix_3 <- expr_matrix
expr_matrix_3 <- t(expr_matrix_2)
view(head(expr_matrix_3))



#---------------------------------------#
                          # Final annotation code of grade normal
norm_data_normal <- rma(data_normal)
expr_matrix <- exprs(norm_data_normal)

# Extract probe-to-gene mapping
probe_ids <- featureNames(norm_data_normal)

head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))
head(probe_ids)

# Extracting data
gene_symbols <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
# gene_entrezid <- mapIds(hta20transcriptcluster.db, keys = probe_ids, column = "ENTREZID", keytype = "PROBEID", multiVals = "first")

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
# gene_entrezid <- gene_entrezid[!is.na(gene_entrezid)]

# Create an annotated expression matrix
expr_matrix <- exprs(norm_data_normal) # Extract expression values

# expr_matrix <- cbind(ENTREZID = gene_entrezid, expr_matrix)
expr_matrix <- cbind(SYMBOL = gene_symbols, expr_matrix)

# convert to data frame and extract duplicates
expr_matrix <- as.data.frame(expr_matrix, stringsAsFactors = FALSE)
duplicates <- expr_matrix %>% filter(duplicated(SYMBOL))
nrow(duplicates)

# Remove duplicates based on SYMBOL
expr_matrix <- expr_matrix %>% distinct(SYMBOL, .keep_all = TRUE)

# Set SYMBOL column as row names
rownames(expr_matrix) <- expr_matrix$SYMBOL  

# Remove the SYMBOL column since it's now in row names
expr_matrix$SYMBOL <- NULL  

# View first few rows
expr_matrix_normal <- expr_matrix
expr_matrix_normal <- t(expr_matrix_normal)
#view(head(expr_matrix_normal))

#-----------------------------------------------------------------------------------#
# adding a column to the datasets. which grade they are. 

expr_matrix_0 <- as.data.frame(expr_matrix_0)
expr_matrix_2 <- as.data.frame(expr_matrix_2)
expr_matrix_3 <- as.data.frame(expr_matrix_3)
expr_matrix_normal <- as.data.frame(expr_matrix_normal)

#-----------------------------------------------------------------------------------#
# Loading all hub genes of all grades. 

hub_genes_0 <- read.csv("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/ppi_interactions/grade_0_hubgenes.csv")
hub_genes_2 <- read.csv("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/ppi_interactions/grade_2_hubgenes.csv")
hub_genes_3 <- read.csv("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/ppi_interactions/grade_3_hubgenes.csv")
view(hub_genes_0)
view(hub_genes_2)
view(hub_genes_3$HubGenes)

hub_genes_0 <- hub_genes_0 %>% select(-Name, -Grade, -Score)
hub_genes_2 <- hub_genes_2 %>% select(-Name,  -Score)
hub_genes_3 <- hub_genes_3 %>% select(-Name,  -Score)


# Find common genes
common_genes_0 <- intersect(hub_genes_0$HubGenes, colnames(expr_matrix_0))
common_genes_0

common_genes_2 <- intersect(hub_genes_2$HubGenes, colnames(expr_matrix_2))
common_genes_2

common_genes_3 <- intersect(hub_genes_3$HubGenes, colnames(expr_matrix_3))
common_genes_3


# Subset the expression data to keep only common genes
expr_matrix_0 <- expr_matrix_0[, colnames(expr_matrix_0) %in% common_genes_0]
expr_matrix_2 <- expr_matrix_2[, colnames(expr_matrix_2) %in% common_genes_2]
expr_matrix_3 <- expr_matrix_3[, colnames(expr_matrix_3) %in% common_genes_3]


expr_matrix_0 <- as.data.frame(expr_matrix_0)
expr_matrix_2 <- as.data.frame(expr_matrix_2)
expr_matrix_3 <- as.data.frame(expr_matrix_3)


expr_matrix_0$Grade <- 0
expr_matrix_2$Grade <- 2
expr_matrix_3$Grade <- 3
expr_matrix_normal$Grade <- 4

view(expr_matrix_0)
view(expr_matrix_2)
view(expr_matrix_3)
view(expr_matrix_normal)

setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs")


write.csv(expr_matrix_0, "hub_genes_0.csv")
# write.csv(expr_matrix_2, "hub_genes_2.csv")
# write.csv(expr_matrix_3, "hub_genes_3.csv")
