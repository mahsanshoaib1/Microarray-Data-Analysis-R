# Pathway Enrichment Analysis
################################################################################

## Details
# we need three things: 1) list of DEGs, 2) background genes(all genes in the organism) 3) list of Gene Sets
# we get list of over represented pathways. 
# There are different websites having gene sets: 1) Ontology. 2) KEGG 3) Reactome

# Before applying PEA we do Over representation Analysis. (ORA)
# ORA is removing DEGs that are not in our threshold. p_value<0.05 and |FC| > 2.
# There are some problems with ORA like what if a gene has p_value=0.049. So, rather than ORA we perform Functional Class Scoring Method(FSC)

# Functional Class Scoring Method(FSC):
# There are different methods to perform FSC): 1) Gene Set Enrichment Analysis (GSEA) 2) Gene Set Variation Analysis. 3) Parametric Analysis of Gene Set Enrichment

## GSEA
# In GSEA, rather then applying threshold to the DEGs, we rank them by level of differential expression. 
# p_value tells how significant the chanage is. 
# The log-2 fold change tells the direction and strength of the changes. Either genes are upregulated or down regulated. 
# The formula to rank values is: ranking = sing(FC)* - log10(p-value)
# it ranks by direction and significance. top genes are significant + upregulated, bottom genes are significant + down regulated, centeral genes are non-significant

##Fisher's exact test
# it is another statistical test used to check if DEGs are enriched in a pathway or not. 
# it also gives a p-value. low p-values means the gene is over represented means involve in pathway. 

##Kolmogorov-Smirnov(kl test)
# The mostly used statistical test use to compare DEGs and Geneset is Kolmogorov-Smirnov(kl test)
# test will give a p value for each pathway. We will correct it for multuple testing.   

## Multiple testing
# then we perform multiple testing to get relevant pathways and genes only. 
# A common multiple testing method is BENJAMINI-HOCHBERG correction. 


##results we get in gene set enrichment analysis
# Enrichment Score(ES): The maximum daviation from zero. explained below in enrichment score plot. 
# Normalized Enrichment Score(NES):
# False Discovery Rate(FDR): 
# nominal p_value:


# Enrichment Score plot/Hiking mountain: it checks each gene. if the gene is present in the geneset, the hike goes up, if it is not present the hike goes down. We want the mountain to be as tall as possible. The maximum hike is th enrichment score or (ES). The are the most important genes in that pathway.
# Leading Edge subest: The leading edge subset of a geneset is the subset of gene members that contribute most of the ES. In plot, they appear in the ranked list before the peak score.  
# if the goal is to compare different pathways we dont use ES, we use NES. 

# for multiple genesets dont use p value, must correct gene set size and multiple hypothesis testing. 

# q value = False Discovery Rate(FDR) are corrected for multiple testing. its probability that a gene set with a certain enrichment score is false positive 
# q value should be samller than 0.25. q_value < 0.25 

# There are different websites having gene sets: 1) Gene Ontology. 2) KEGG 3) Reactome

# Gene Ontology: Focuses on biological processes. 
# KEGG: Focuses on metabolic pathways. 
# Reactome: Human Molecular Pathways. 

# Gene Sets are not restricted to function. There are also gene sets of diseases, gene tissues, transciption transcriptor factors. 

####################################################################################################################
# method/video - 1
# loading libraries
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggupset)
library(pheatmap)
#library(topGO)
library(hgu133a.db)  # Human gene annotation
library(DOSE)           # For enrichment analysis
library(enrichplot)     # For visualization
library(ggplot2)
library(readxl)

# setting paths
setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/pathway")
path <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/pathway" # input path, where your data is located


#------------------------------------------------------------------------#


# loading Gene Data
# Load the CSV file containing common genes (ensure the file is in the working directory)
gene_data <- read_excel("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/overlapped/Grade_2_Common_genes.xlsx")
#view(gene_data)
# removing duplicate genes
gene_data <- gene_data[!duplicated(gene_data$GENENAME),] 

# extracting logFC column out of all the data
gene_list <- gene_data$logFC
#view(gene_list)

# genes <- as.data.frame(logFC) # converting logFC data into dataframe

# making GENENAME as rowname for logFC dataframe
names(gene_list) <- gene_data$GENENAME 
#view(gene_list)

# sorting the dataframes
gene_list <- sort.DataFrame(gene_list, decreasing = TRUE)
#view(gene_list)

# setting organism
library(org.Hs.eg.db, character.only == TRUE)

# specifying gene ontology code term
# BP = Biological processes
# CC = Cell Cycle
# MF = Molecular Function
# GSEA = Gene Set Enrichment analysis

gse_go <- gseGO(
  geneList = gene_list, 
  ont = 'All',
  keyType = 'SYMBOL', 
  minGSSize = 3, 
  maxGSSize = 100, 
  pvalueCutoff = 0.05,
  verbose = TRUE, 
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = 'none'
)

view(gse_go)

write.csv(gse_go, "Grade_2_gse_go_results.csv")


# summary <- summary(gse_go)
# view(summary)

#-----------------------------#
# KEGG
# converting gene symbol to entrezid because kegg takes ENTREZID onlt
# Convert SYMBOLs to Entrez IDs
converted_genes <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(converted_genes)
# head(entres)

# Merge with original logFC values
converted_genes$logFC <- gene_list[converted_genes$SYMBOL]
view(converted_genes)

# Create a named vector for gseKEGG (Entrez IDs as names)
gene_list_entrez <- setNames(converted_genes$logFC, converted_genes$ENTREZID)

# Ensure it's sorted in decreasing order
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
head(gene_list_entrez)
head(gene_list)

gse_kegg <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  exponent = 1,
  minGSSize = 3, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05,
  verbose = TRUE, 
  pAdjustMethod = 'BH',
)

# view(gse_kegg)
write.csv(gse_kegg, "Grade_2_gse_kegg_results.csv")

# summary <- summary(gse_kegg)
# view(summary)

# enrichKEGG(
#   gene_list_entrez,
#   organism = "hsa",
#   keyType = "ncbi-geneid",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.2,
#   use_internal_data = FALSE
# )
# 


#----------------------------------------------------------------------------#

####################################################################################################################
# method/video - 2






