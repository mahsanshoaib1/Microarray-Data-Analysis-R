# import pandas as pd
# from matplotlib_venn import venn2
# import matplotlib.pyplot as plt

# # Load both Excel files
# alz_df = pd.read_excel("C:\\Users\\urbaz\\OneDrive\\Desktop\\R\\GEOD-1297\\alzheimer_data.xlsx")
# pd_df = pd.read_excel("C:\\Users\\urbaz\\OneDrive\\Desktop\\R\\E-GEOD-20168\\Parkinson_data.xlsx")

# # Make sure to use the correct column name containing gene symbols
# # Example: if the gene column is "Gene Symbol", change it here.
# # alz_genes = set(alz_df['Gene Symbol'])  # Replace 'Gene Symbol' with actual column name if different
# # pd_genes = set(pd_df['Gene Symbol'])

# # Plot Venn Diagram
# plt.figure(figsize=(8, 8))
# venn = venn2([alz_genes, pd_genes], ('Alzheimer\'s DEGs', 'Parkinson\'s DEGs'))

# # Add title (optional)
# plt.title('Overlap Between Alzheimer\'s and Parkinson\'s DEGs')

# # Show plot
# plt.show()

# # Optional: If you want to also save the plot
# # plt.savefig('DEGs_Overlap_Venn.png', dpi=300)



import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Load TNBC GRADE 2 DEGs file
grade_2_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\grade_2\\updated_Grade_2_TNBC_DEG_results.csv')
grade_2_genes = set(grade_2_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for TNBC grade 2
grade_2_gene_dict = dict(zip(grade_2_df['SYMBOL'], grade_2_df['GENENAME']))


# Load Ovarian DEGs file
ovarian_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\ovarian\\updated_ovarian_DEG_results.csv')
ovarian_genes = set(ovarian_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for Parkinson’s
ovarian_gene_dict = dict(zip(ovarian_df['SYMBOL'], ovarian_df['GENENAME']))


# Find common genes
common_genes_2 = grade_2_genes.intersection(ovarian_genes)


# Create a list to store common genes with their full names
common_gene_list_2 = []

for gene in common_genes_2:
    genename = grade_2_gene_dict.get(gene, ovarian_gene_dict.get(gene, 'Unknown'))  # Prefer disease's name if available
    common_gene_list_2.append((gene, genename))


# Print or save the common genes list
print(type(common_gene_list_2))
print(f"Number of common genes in TNBC Grade 2 and ovarian cancer: {len(common_gene_list_2)}")
for symbol, genename in common_gene_list_2:
    print(f"{symbol}\t{genename}")

# adding the other coloumns
common_df = pd.DataFrame(common_gene_list_2, columns=['SYMBOL', 'GENENAME'])

# Plot Venn diagram
plt.figure(figsize=(8, 6))
venn2([grade_2_genes, ovarian_genes], ('TNBC Grade 2 DEGs', 'ovarian DEGs'))
plt.title('Overlap of TNBC Grade 2 DEGs between TNBC and ovarian')
plt.show()
plt.savefig('Grade_2_DEGs_Overlap_Venn.png', dpi=300)


# Extract required columns from grade_2_df
selected_columns = ['SYMBOL', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
grade_2_selected = grade_2_df[selected_columns]

# Merge selected columns into common_df based on SYMBOL
common_df = common_df.merge(grade_2_selected, on='SYMBOL', how='left')

# Optional: Save to Excel
common_df.to_excel("D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\overlapped\\Grade_2_Common_genes.xlsx", index=False)





################################################################################
################################################################################


# for grade 3
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt


# Load TNBC GRADE 3 DEGs file
grade_3_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\grade_3\\updated_Grade_3_TNBC_DEG_results.csv')
grade_3_genes = set(grade_3_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for TNBC grade 3
grade_3_gene_dict = dict(zip(grade_3_df['SYMBOL'], grade_3_df['GENENAME']))


# Load Ovarian DEGs file
ovarian_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\ovarian\\updated_ovarian_DEG_results.csv')
ovarian_genes = set(ovarian_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for Parkinson’s
ovarian_gene_dict = dict(zip(ovarian_df['SYMBOL'], ovarian_df['GENENAME']))


# Find common genes
common_genes_3 = grade_3_genes.intersection(ovarian_genes)


# Create a list to store common genes with their full names
common_gene_list_3 = []

for gene in common_genes_3:
    genename = grade_3_gene_dict.get(gene, ovarian_gene_dict.get(gene, 'Unknown'))  # Prefer Alzheimer's name if available
    common_gene_list_3.append((gene, genename))


# Print or save the common genes list
print(f"Number of common genes in TNBC Grade 3 and ovarian cancer: {len(common_gene_list_3)}")
for symbol, genename in common_gene_list_3:
    print(f"{symbol}\t{genename}")

# adding the other coloumns
common_df = pd.DataFrame(common_gene_list_3, columns=['SYMBOL', 'GENENAME'])

# Plot Venn diagram
plt.figure(figsize=(8, 6))
venn2([grade_3_genes, ovarian_genes], ('Grade 3 TNBC DEGs', 'ovarian DEGs'))
plt.title('Overlap of DEGs between TNBC Grade 3 and ovarian')
plt.show()
plt.savefig('Grade_3_DEGs_Overlap_Venn.png', dpi=300)


# Extract required columns from grade_3_df
selected_columns = ['SYMBOL', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
grade_3_selected = grade_3_df[selected_columns]

# Merge selected columns into common_df based on SYMBOL
common_df = common_df.merge(grade_3_selected, on='SYMBOL', how='left')

# Optional: Save to Excel
common_df.to_excel("D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\overlapped\\Grade_3_Common_genes.xlsx", index=False)





################################################################################
################################################################################

# for grade 0
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt


# Load TNBC GRADE 3 DEGs file
grade_0_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\grade_0\\updated_Grade_0_TNBC_DEG_results.csv')
grade_0_genes = set(grade_0_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for TNBC grade 0
grade_0_gene_dict = dict(zip(grade_0_df['SYMBOL'], grade_0_df['GENENAME']))


# Load Ovarian DEGs file
ovarian_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\ovarian\\updated_ovarian_DEG_results.csv')
ovarian_genes = set(ovarian_df['SYMBOL'])
# Create a dictionary mapping SYMBOL to GENENAME for Parkinson’s
ovarian_gene_dict = dict(zip(ovarian_df['SYMBOL'], ovarian_df['GENENAME']))


# Find common genes
common_genes_0 = grade_0_genes.intersection(ovarian_genes)



# Create a list to store common genes with their full names
common_gene_list_0 = []

for gene in common_genes_0:
    genename = grade_0_gene_dict.get(gene, ovarian_gene_dict.get(gene, 'Unknown'))  # Prefer Alzheimer's name if available
    common_gene_list_0.append((gene, genename))


# Print or save the common genes list
print(f"Number of common genes in TNBC Grade 0 and ovarian: {len(common_gene_list_0)}")
for symbol, genename in common_gene_list_0:
    print(f"{symbol}\t{genename}")


# adding the other coloumns
common_df = pd.DataFrame(common_gene_list_0, columns=['SYMBOL', 'GENENAME'])

# Plot Venn diagram
plt.figure(figsize=(8, 6))
venn2([grade_0_genes, ovarian_genes], ('Grade 0 TNBC DEGs', 'ovarian DEGs'))
plt.title('Overlap of DEGs between Grade 0 TNBC and ovarian')
plt.show()
plt.savefig('Grade_0_DEGs_Overlap_Venn.png', dpi=300)

import logging

# Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create a file handler
file_handler = logging.FileHandler('app.log')
file_handler.setLevel(logging.INFO)

# Create a console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Create a formatter and set it for the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# ... rest of the code ...

# Load TNBC GRADE 2 DEGs file
logger.info('Loading TNBC GRADE 2 DEGs file')
grade_2_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\grade_2\\updated_Grade_2_TNBC_DEG_results.csv')
logger.info('TNBC GRADE 2 DEGs file loaded')

# ... rest of the code ...

# Load Ovarian DEGs file
logger.info('Loading Ovarian DEGs file')
ovarian_df = pd.read_csv('D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\ovarian\\updated_ovarian_DEG_results.csv')
logger.info('Ovarian DEGs file loaded')

# ... rest of the code ...

# Find common genes
logger.info('Finding common genes')
common_genes_2 = grade_2_genes.intersection(ovarian_genes)
logger.info('Common genes found')

# ... rest of the code ...

# Plot Venn diagram
logger.info('Plotting Venn diagram')
plt.figure(figsize=(8, 6))
venn2([grade_2_genes, ovarian_genes], ('TNBC Grade 2 DEGs', 'ovarian DEGs'))
plt.title('Overlap of TNBC Grade 2 DEGs between TNBC and ovarian')
plt.show()
plt.savefig('Grade_2_DEGs_Overlap_Venn.png', dpi=300)
logger.info('Venn diagram plotted')

# ... rest of the code ...

# Extract required columns from grade_2_df
logger.info('Extracting required columns from grade_2_df')
selected_columns = ['SYMBOL', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
grade_2_selected = grade_2_df[selected_columns]
logger.info('Required columns extracted')

# ... rest of the code ...

# Merge selected columns into common_df based on SYMBOL
logger.info('Merging selected columns into common_df')
common_df = common_df.merge(grade_2_selected, on='SYMBOL', how='left')
logger.info('Selected columns merged')

# ... rest of the code ...

# Optional: Save to Excel
logger.info('Saving to Excel')
common_df.to_excel("D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\overlapped\\Grade_2_Common_genes.xlsx", index=False)
logger.info('Excel file saved')
# Extract required columns from grade_0_df
selected_columns = ['SYMBOL', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
grade_0_selected = grade_0_df[selected_columns]

# Merge selected columns into common_df based on SYMBOL
common_df = common_df.merge(grade_0_selected, on='SYMBOL', how='left')

# Optional: Save to Excel
common_df.to_excel("D:\\study\\lab_projects\\Dr_Zubair\\Breast _cancer\\practice_practical\\Selected Datasets\\outputs\\overlapped\\Grade_0_Common_genes.xlsx", index=False)




