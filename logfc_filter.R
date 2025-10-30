# Install required packages if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Load required libraries
library(ggplot2)
library(pheatmap)

# Placeholder: Load expression data and sample annotations
# REPLACE THIS WITH YOUR ACTUAL EXPRESSION DATA
# Assume expr_data is a matrix/data frame with rows as genes (matching GeneID) and columns as samples
# Example: expr_data <- read.csv("path/to/your/expression_data.csv", row.names = 1)
# For now, this is a placeholder (you must provide the actual data)
expr_data <- matrix(NA, nrow = 1, ncol = 1) # Placeholder

# Placeholder: Define sample groups for heatmap annotations
# REPLACE THIS WITH YOUR ACTUAL SAMPLE METADATA
# Example: sample_groups <- data.frame(Sample = colnames(expr_data), Condition = c("TNBC", "Normal", ...))
# For now, this is a placeholder
sample_groups <- data.frame(Sample = character(), Condition = character()) # Placeholder

# to apply log fc filter on the differentially expressed genes and generate visualizations

############################### GRADE 0  ###############################

# Set the input and output paths
input_file <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/grade_0/updated_Grade_0_TNBC_DEG_results.csv"
output_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/logfc_output"
output_file <- file.path(output_dir, "filtered_Grade_0_TNBC_DEG_results.csv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the data
deg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Apply filters: |logFC| > 1 and adj.P.Val < 0.05
filtered_degs <- deg_data[abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, ]

# Save the filtered results to the output directory
write.csv(filtered_degs, output_file, row.names = FALSE)

# Print a message to confirm
cat("Filtered DEGs saved to:", output_file, "\n")
cat("Number of DEGs after filtering:", nrow(filtered_degs), "\n")

# Visualization: Volcano Plot for Grade 0
deg_data$Significance <- ifelse(abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, "Significant", "Not Significant")
volcano_plot <- ggplot(deg_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Grade 0 TNBC vs Normal", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
ggsave(file.path(output_dir, "volcano_Grade_0_TNBC.pdf"), volcano_plot, width = 8, height = 6)

# Visualization: Heatmap for Grade 0
if (nrow(filtered_degs) > 0 && exists("expr_data") && nrow(expr_data) > 0) {
  # Subset expression data to filtered DEGs (assuming GeneID column exists)
  deg_expr <- expr_data[rownames(expr_data) %in% filtered_degs$GeneID, , drop = FALSE]
  if (nrow(deg_expr) > 0) {
    # Prepare annotation for samples
    annotation_col <- data.frame(Condition = sample_groups$Condition, row.names = sample_groups$Sample)
    # Generate heatmap
    pdf(file.path(output_dir, "heatmap_Grade_0_TNBC.pdf"), width = 10, height = 8)
    pheatmap(deg_expr, 
             scale = "row", 
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_col = annotation_col,
             main = "Heatmap: Grade 0 TNBC DEGs",
             show_rownames = nrow(deg_expr) <= 50)
    dev.off()
  }
}

############################### GRADE 2  ###############################

# Set the input and output paths
input_file <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/grade_2/updated_Grade_2_TNBC_DEG_results.csv"
output_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/logfc_output"
output_file <- file.path(output_dir, "filtered_Grade_2_TNBC_DEG_results.csv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the data
deg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Apply filters: |logFC| > 1 and adj.P.Val < 0.05
filtered_degs <- deg_data[abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, ]

# Save the filtered results to the output directory
write.csv(filtered_degs, output_file, row.names = FALSE)

# Print a message to confirm
cat("Filtered DEGs saved to:", output_file, "\n")
cat("Number of DEGs after filtering:", nrow(filtered_degs), "\n")

# Visualization: Volcano Plot for Grade 2
deg_data$Significance <- ifelse(abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, "Significant", "Not Significant")
volcano_plot <- ggplot(deg_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Grade 2 TNBC vs Normal", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
ggsave(file.path(output_dir, "volcano_Grade_2_TNBC.pdf"), volcano_plot, width = 8, height = 6)

# Visualization: Heatmap for Grade 2
if (nrow(filtered_degs) > 0 && exists("expr_data") && nrow(expr_data) > 0) {
  deg_expr <- expr_data[rownames(expr_data) %in% filtered_degs$GeneID, , drop = FALSE]
  if (nrow(deg_expr) > 0) {
    annotation_col <- data.frame(Condition = sample_groups$Condition, row.names = sample_groups$Sample)
    pdf(file.path(output_dir, "heatmap_Grade_2_TNBC.pdf"), width = 10, height = 8)
    pheatmap(deg_expr, 
             scale = "row", 
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_col = annotation_col,
             main = "Heatmap: Grade 2 TNBC DEGs",
             show_rownames = nrow(deg_expr) <= 50)
    dev.off()
  }
}

############################### GRADE 3  ###############################

# Set the input and output paths
input_file <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/grade_3/updated_Grade_3_TNBC_DEG_results.csv"
output_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/logfc_output"
output_file <- file.path(output_dir, "filtered_Grade_3_TNBC_DEG_results.csv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the data
deg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Apply filters: |logFC| > 1 and adj.P.Val < 0.05
filtered_degs <- deg_data[abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, ]

# Save the filtered results to the output directory
write.csv(filtered_degs, output_file, row.names = FALSE)

# Print a message to confirm
cat("Filtered DEGs saved to:", output_file, "\n")
cat("Number of DEGs after filtering:", nrow(filtered_degs), "\n")

# Visualization: Volcano Plot for Grade 3
deg_data$Significance <- ifelse(abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, "Significant", "Not Significant")
volcano_plot <- ggplot(deg_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Grade 3 TNBC vs Normal", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
ggsave(file.path(output_dir, "volcano_Grade_3_TNBC.pdf"), volcano_plot, width = 8, height = 6)

# Visualization: Heatmap for Grade 3
if (nrow(filtered_degs) > 0 && exists("expr_data") && nrow(expr_data) > 0) {
  deg_expr <- expr_data[rownames(expr_data) %in% filtered_degs$GeneID, , drop = FALSE]
  if (nrow(deg_expr) > 0) {
    annotation_col <- data.frame(Condition = sample_groups$Condition, row.names = sample_groups$Sample)
    pdf(file.path(output_dir, "heatmap_Grade_3_TNBC.pdf"), width = 10, height = 8)
    pheatmap(deg_expr, 
             scale = "row", 
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_col = annotation_col,
             main = "Heatmap: Grade 3 TNBC DEGs",
             show_rownames = nrow(deg_expr) <= 50)
    dev.off()
  }
}

############################### GRADE ovarian  ###############################

# Set the input and output paths
input_file <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/ovarian/updated_ovarian_DEG_results.csv"
output_dir <- "D:/study/lab_projects/Dr_Zubair/Breast _cancer/practice_practical/Selected Datasets/outputs/logfc_output"
output_file <- file.path(output_dir, "filtered_ovarian_TNBC_DEG_results.csv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the data
deg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Apply filters: |logFC| > 1 and adj.P.Val < 0.05
filtered_degs <- deg_data[abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, ]

# Save the filtered results to the output directory
write.csv(filtered_degs, output_file, row.names = FALSE)

# Print a message to confirm
cat("Filtered DEGs saved to:", output_file, "\n")
cat("Number of DEGs after filtering:", nrow(filtered_degs), "\n")

# Visualization: Volcano Plot for Ovarian
deg_data$Significance <- ifelse(abs(deg_data$logFC) > 1 & deg_data$adj.P.Val < 0.05, "Significant", "Not Significant")
volcano_plot <- ggplot(deg_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Ovarian vs Control", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
ggsave(file.path(output_dir, "volcano_ovarian.pdf"), volcano_plot, width = 8, height = 6)

# Visualization: Heatmap for Ovarian
if (nrow(filtered_degs) > 0 && exists("expr_data") && nrow(expr_data) > 0) {
  deg_expr <- expr_data[rownames(expr_data) %in% filtered_degs$GeneID, , drop = FALSE]
  if (nrow(deg_expr) > 0) {
    annotation_col <- data.frame(Condition = sample_groups$Condition, row.names = sample_groups$Sample)
    pdf(file.path(output_dir, "heatmap_ovarian.pdf"), width = 10, height = 8)
    pheatmap(deg_expr, 
             scale = "row", 
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_col = annotation_col,
             main = "Heatmap: Ovarian DEGs",
             show_rownames = nrow(deg_expr) <= 50)
    dev.off()
  }
}