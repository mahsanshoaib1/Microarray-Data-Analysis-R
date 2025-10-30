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
setwd("D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/overlapped")
input_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/overlapped"
output_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/pathway_output"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a log file for KEGG diagnostics
log_file <- file.path(output_dir, "kegg_diagnostic_log.txt")
cat("KEGG Analysis Diagnostic Log\n", file = log_file)

# Define input files and corresponding output file names
datasets <- list(
  list(
    input_file = "Grade_0_Common_genes.xlsx",
    go_output = "Grade_0_gse_go_results.csv",
    kegg_output = "Grade_0_gse_kegg_results.csv"
  ),
  list(
    input_file = "Grade_2_Common_genes.xlsx",
    go_output = "Grade_2_gse_go_results.csv",
    kegg_output = "Grade_2_gse_kegg_results.csv"
  ),
  list(
    input_file = "Grade_3_Common_genes.xlsx",
    go_output = "Grade_3_gse_go_results.csv",
    kegg_output = "Grade_3_gse_kegg_results.csv"
  )
)

# Loop through each dataset
for (dataset in datasets) {
  # Load the Excel file containing common genes
  gene_data <- readxl::read_excel(file.path(input_dir, dataset$input_file))
  
  # Check if SYMBOL column exists; if not, assume GENENAME is SYMBOL
  if (!"SYMBOL" %in% colnames(gene_data)) {
    if ("GENENAME" %in% colnames(gene_data)) {
      gene_data$SYMBOL <- gene_data$GENENAME
      cat("Using GENENAME as SYMBOL for dataset:", dataset$input_file, "\n")
    } else {
      stop("Neither SYMBOL nor GENENAME found in ", dataset$input_file)
    }
  }
  
  # Remove duplicate genes based on SYMBOL
  gene_data <- gene_data[!duplicated(gene_data$SYMBOL), ]
  
  # Extract logFC column and name it with gene symbols
  if (!"logFC" %in% colnames(gene_data)) {
    stop("logFC column not found in ", dataset$input_file)
  }
  gene_list <- gene_data$logFC
  names(gene_list) <- gene_data$SYMBOL
  
  # Sort the gene list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # GSEA for Gene Ontology (GO)
  gse_go <- tryCatch({
    gseGO(
      geneList = gene_list, 
      ont = 'ALL',
      keyType = 'SYMBOL', 
      minGSSize = 3, 
      maxGSSize = 100, 
      pvalueCutoff = 0.05,
      verbose = TRUE, 
      OrgDb = org.Hs.eg.db,
      pAdjustMethod = 'none'
    )
  }, error = function(e) {
    cat("Error in gseGO for dataset:", dataset$input_file, "\n", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Save GO results if successful
  if (!is.null(gse_go)) {
    write.csv(gse_go, file.path(output_dir, dataset$go_output), row.names = FALSE)
    cat("GO results saved to:", file.path(output_dir, dataset$go_output), "\n")
  } else {
    cat("No GO results saved for:", dataset$input_file, "\n")
  }
  
  # KEGG Analysis (Updated)
  # Log dataset details
  cat("Processing KEGG for dataset:", dataset$input_file, "\n", file = log_file, append = TRUE)
  cat("Number of genes loaded:", nrow(gene_data), "\n", file = log_file, append = TRUE)
  cat("Columns in dataset:", paste(colnames(gene_data), collapse = ", "), "\n", file = log_file, append = TRUE)
  cat("Number of genes after deduplication:", nrow(gene_data), "\n", file = log_file, append = TRUE)
  
  # Handle NA logFC values
  na_count <- sum(is.na(gene_list))
  cat("Number of NA logFC values:", na_count, "\n", file = log_file, append = TRUE)
  gene_list <- gene_list[!is.na(gene_list)]
  cat("Number of genes in gene_list:", length(gene_list), "\n", file = log_file, append = TRUE)
  
  if (length(gene_list) == 0) {
    cat("Empty gene list after removing NA logFC\n", file = log_file, append = TRUE)
    cat("Empty gene list for KEGG in ", dataset$input_file, "\n")
    next
  }
  
  # Convert SYMBOLs to Entrez IDs for KEGG
  converted_genes <- tryCatch({
    bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat("Error in bitr:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    cat("Error in bitr for dataset:", dataset$input_file, "\n", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(converted_genes)) {
    cat("No genes mapped to Entrez IDs\n", file = log_file, append = TRUE)
    cat("No KEGG analysis performed for:", dataset$input_file, "due to gene mapping failure\n")
    next
  }
  
  # Log mapping results
  cat("Number of genes mapped to Entrez IDs:", nrow(converted_genes), "\n", file = log_file, append = TRUE)
  
  # Merge with original logFC values
  converted_genes$logFC <- gene_list[converted_genes$SYMBOL]
  cat("Number of genes with logFC after mapping:", sum(!is.na(converted_genes$logFC)), "\n", file = log_file, append = TRUE)
  
  # Create a named vector for gseKEGG (Entrez IDs as names)
  gene_list_entrez <- setNames(converted_genes$logFC, converted_genes$ENTREZID)
  gene_list_entrez <- gene_list_entrez[!is.na(gene_list_entrez)]
  
  # Sort in decreasing order
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
  cat("Number of genes in gene_list_entrez:", length(gene_list_entrez), "\n", file = log_file, append = TRUE)
  
  # Skip if too few genes
  if (length(gene_list_entrez) < 1) {
    cat("Too few genes for KEGG analysis (< 1)\n", file = log_file, append = TRUE)
    cat("Too few genes for KEGG analysis in ", dataset$input_file, "\n")
    next
  }
  
  # GSEA for KEGG
  gse_kegg <- tryCatch({
    gseKEGG(
      geneList = gene_list_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      exponent = 1,
      minGSSize = 1,  # Allow small gene sets
      maxGSSize = 1000,  # Include larger pathways
      pvalueCutoff = 1,  # Lenient to capture all pathways
      verbose = TRUE,
      pAdjustMethod = 'BH',
      nPermSimple = 1000  # Disable parallel processing
    )
  }, error = function(e) {
    cat("Error in gseKEGG:", conditionMessage(e), "\n", file = log_file, append = TRUE)
    cat("Error in gseKEGG for dataset:", dataset$input_file, "\n", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Save KEGG results if successful
  if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
    write.csv(gse_kegg, file.path(output_dir, dataset$kegg_output), row.names = FALSE)
    cat("KEGG results saved to:", dataset$kegg_output, "\n", file = log_file, append = TRUE)
    cat("KEGG results saved to:", file.path(output_dir, dataset$kegg_output), "\n")
  } else {
    cat("No KEGG results saved (empty or null)\n", file = log_file, append = TRUE)
    cat("No KEGG results saved for:", dataset$input_file, "\n")
  }
  
  # Print confirmation
  cat("Processed dataset:", dataset$input_file, "\n\n")
}