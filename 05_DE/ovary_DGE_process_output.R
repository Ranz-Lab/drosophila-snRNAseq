from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import glob
import os

# Define the predetermined order of cell types
cell_type_order = [
    "Germarium 1 and 2a","Germarium 2a 2b","Germarium 2b 3","Stalk & Polar Cells","Follicle Stem Cells & preFCs","Early Follicle Cells",
    "Mitotic Follicle Cells Stage 1 to 5","Post Mitotic Follicle Cells Stage 6","Vitellogenic MBFCs Stage 7",
    "Vitellogenic MBFCs Stage 8","Vitellogenic MBFCs Stage 9 10A","Choriogenic MBFCs Stage 12","Choriogenic MBFCs Stage 14",
    "Terminal Corpus Luteum Cells","Stretch Cells","Ovarian Sheath Muscle","Oviduct"
]

# Define the order of comparisons
comparison_order = [
    "A4 vs w501",
    "ISO1 vs w501",
    "ISO1 vs A4"
]

# Step 1: Read and combine input files into a single DataFrame
file_pattern = "DGE_*_ovary_*.p-val.oct8.txt"
file_list = glob.glob(file_pattern)
combined_data = []
gene_to_id = {}

print("Reading input files...")
for file_name in file_list:
    parts = file_name.split('_')
    comparison = parts[1].replace('vs', ' vs ')
    celltype = ' '.join(parts[3:]).split('.p')[0].replace('-', ' ')
    
    print(f"Processing file: {file_name}")
    df = pd.read_csv(file_name, sep='\t')
    
    if comparison not in comparison_order:
        continue
    
    df['comparison'] = comparison
    df['celltype'] = celltype
    
    for gene in df['gene']:
        if gene not in gene_to_id:
            gene_to_id[gene] = len(gene_to_id) + 1
    df['gene_id'] = df['gene'].map(gene_to_id)
    
    combined_data.append(df)

combined_df = pd.concat(combined_data)
print(f"Combined DataFrame shape: {combined_df.shape}")
print("Combined DataFrame head:\n", combined_df.head())

# Step 2: Reorder the DataFrame based on the predetermined order
combined_df['celltype'] = pd.Categorical(combined_df['celltype'], categories=cell_type_order, ordered=True)
combined_df['comparison'] = pd.Categorical(combined_df['comparison'], categories=comparison_order, ordered=True)
combined_df = combined_df.sort_values(by=['comparison', 'celltype'])


# Step 3: Calculate adjusted p-values using Benjamini-Hochberg for Expressed genes only
print("Calculating adjusted p-values for Expressed genes...")

# Filter for rows where 'Expressed' column is 'E'
expressed_mask = combined_df['Expressed'] == 'E'
expressed_genes_df = combined_df[expressed_mask]

# Perform the B/H correction on all expressed genes
all_p_values = expressed_genes_df['p_val_adj'].values
rejected, adjusted_p_values, _, _ = multipletests(all_p_values, method='fdr_bh')

# Create a new column for adjusted p-values, initially filled with '1.0' for all NE genes
combined_df['adj_p_val'] = 1.0

# Update the adj_p_val column for expressed genes
combined_df.loc[expressed_mask, 'adj_p_val'] = adjusted_p_values.astype(float)

# see the adjusted p-values in scientific notation:
combined_df['adj_p_val'] = combined_df['adj_p_val'].map(lambda x: f"{x:.6e}")  # Adjust the number of decimal places as needed

print("Adjusted p-values added to DataFrame")
print("\nSample of expressed genes:")
print(combined_df[combined_df['Expressed'] == 'E'].head(15))
print("\nSample of non-expressed genes:")
print(combined_df[combined_df['Expressed'] == 'NE'].head(15))

# Step 4: Define genes as "DE" or "nDE", with "NA" for non-expressed genes
def define_de_status(row):
    if row['Expressed'] != 'E':
        return np.nan  # Not expressed, set to NA
    elif row['adj_p_val_float'] < 0.01 and abs(row['avg_log2FC']) > 1:
        return 1 if row['avg_log2FC'] > 0 else -1  # 1 for upregulation, -1 for downregulation
    else:
        return 0  # nDE (Not Differentially Expressed)

# Convert adj_p_val back to float for comparison, handling 'X' and 'nan' values
combined_df['adj_p_val_float'] = pd.to_numeric(combined_df['adj_p_val'].replace('X', np.nan), errors='coerce')

print("Defining DE status, with NA for non-expressed genes...")
combined_df['DE_status'] = combined_df.apply(define_de_status, axis=1)

# Convert the DE_status column to Int64 type to allow integers and NaNs
combined_df['DE_status'] = combined_df['DE_status'].astype('Int64')

print("\nSample of expressed genes:")
print(combined_df[combined_df['Expressed'] == 'E'].shape)
print(combined_df[combined_df['Expressed'] == 'E'].head(15))
print("\nSample of non-expressed genes:")
print(combined_df[combined_df['Expressed'] == 'NE'].shape)
print(combined_df[combined_df['Expressed'] == 'NE'].head(15))

# Step 5: Compare with original list, but only for expressed genes
def original_de_status(regulation):
    if regulation == 'Up':
        return 1
    elif regulation == 'Down':
        return -1
    else:
        return 0

print("Performing comparison with original DE status for expressed genes only...")

# Initialize the new columns with placeholders
combined_df['original_DE_status'] = 'X'
combined_df['comparison_result'] = 'X'

# Filter for expressed genes
expressed_mask = combined_df['Expressed'] == 'E'

# Apply the original DE status function and comparison only for expressed genes
combined_df.loc[expressed_mask, 'original_DE_status'] = combined_df.loc[expressed_mask, 'regulation'].apply(original_de_status)

# Compare the calculated DE_status with the original_DE_status for expressed genes
combined_df.loc[expressed_mask, 'comparison_result'] = np.where(
    combined_df.loc[expressed_mask, 'DE_status'] == combined_df.loc[expressed_mask, 'original_DE_status'], 
    "OK", 
    "discrepancy"
)

print("Comparison result added to DataFrame")
print("\nSample of expressed genes:")
print(combined_df[combined_df['Expressed'] == 'E'].tail(15))
print("\nSample of non-expressed genes:")
print(combined_df[combined_df['Expressed'] == 'NE'].head(15))

# Step 6: Create comparison summary, adding 'X' for NE genes and retaining current summary for E genes

print("Creating comparison summary with 'X' for NE genes...")

# Create a DataFrame for DE binary values with conditions for 'Expressed' column
de_binary_df = pd.DataFrame({
    'gene': combined_df['gene'],
    'comparison': combined_df['comparison'],
    'celltype': combined_df['celltype'],
    'DE_binary': np.where(
        combined_df['Expressed'] == 'E', 
        combined_df['DE_status'].astype(str),  # Convert to string for E genes
        'X'  # Use 'X' for NE genes
    )
})

print("DE Binary DataFrame head:\n", de_binary_df.head())

# Pivot to create comparison summary
comparison_summary = de_binary_df.pivot_table(
    index='gene',
    columns=['comparison', 'celltype'],
    values='DE_binary',
    aggfunc='first',
    fill_value='0'
)

print("Comparison summary before reindexing:\n", comparison_summary.head())

# Reindex to ensure correct order
comparison_summary = comparison_summary.reindex(columns=pd.MultiIndex.from_product([comparison_order, cell_type_order], names=['comparison', 'celltype']), fill_value='0')

# Simple string replacement to swap out -1.0, 1.0, and 0.0 to -1, 1, and 0
#comparison_summary = comparison_summary.replace({'-1.0': '-1', '1.0': '1', '0.0': '0'})

print("Comparison summary shape after reindexing:", comparison_summary.shape)
print("Comparison summary head after reindexing:\n", comparison_summary.tail(25))

# Step 7: Determine combination of DEGs and calculate percentages
print("Determining combinations of DEGs...")
comparison_summary['combo_A4_vs_w501'] = comparison_summary[comparison_order[0]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
comparison_summary['combo_ISO1_vs_w501'] = comparison_summary[comparison_order[1]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
comparison_summary['combo_ISO1_vs_A4'] = comparison_summary[comparison_order[2]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

combo_counts_A4_vs_w501 = comparison_summary['combo_A4_vs_w501'].value_counts(normalize=True) * 100
combo_counts_ISO1_vs_w501 = comparison_summary['combo_ISO1_vs_w501'].value_counts(normalize=True) * 100
combo_counts_ISO1_vs_A4 = comparison_summary['combo_ISO1_vs_A4'].value_counts(normalize=True) * 100

print("Combo counts for A4 vs w501:\n", combo_counts_A4_vs_w501)
print("Combo counts for ISO1 vs w501:\n", combo_counts_ISO1_vs_w501)
print("Combo counts for ISO1 vs A4:\n", combo_counts_ISO1_vs_A4)


# Step 8: Create final combinations
print("Creating final combinations...")
# Convert 'DE_status' to string type, then fill NA values with 'X'
combined_df['DE_status'] = combined_df['DE_status'].astype('string').fillna('X')

# Proceed with the groupby operation
final_combinations = combined_df.groupby('gene')['DE_status'].apply(
    lambda x: ''.join(['1' if status == '1' else '-1' if status == '-1' else '0' if status == '0' else 'X' for status in x])
)

final_combinations.columns = ['gene', 'DE_combination']
print("Final combinations head:\n", final_combinations.head(25))


# Step 9: Consider a gene as '1' if DE in at least one cell type or as 'X' if all are 'X'
print("Considering DE in at least one cell type...")

at_least_one_de = combined_df.groupby('gene')['DE_status'].apply(
    lambda x: 'X' if all(x == 'X') else ('DE' if any(pd.to_numeric(x, errors='coerce').fillna(0) != 0) else 'nDE')
).reset_index()

# Set column names
at_least_one_de.columns = ['gene', 'at_least_one_DE']

# Convert 'final_combinations' from a Series to a DataFrame if necessary
final_combinations = final_combinations.reset_index()
final_combinations.columns = ['gene', 'DE_combination']

# Now, merge with 'at_least_one_de'
final_combinations = final_combinations.merge(at_least_one_de, on='gene', how='left')

print("Final combinations with at least one DE:\n", final_combinations.head(25))


# Save results to tab-separated text files
output_dir = os.getcwd()
print("Saving results to files...")
combined_df.to_csv(os.path.join(output_dir, "combined_results.txt"), sep='\t', index=False)
comparison_summary.to_csv(os.path.join(output_dir, "comparison_summary.txt"), sep='\t')
final_combinations.to_csv(os.path.join(output_dir, "final_combinations.txt"), sep='\t', index=False)
combo_counts_A4_vs_w501.to_csv(os.path.join(output_dir, "combo_counts_A4_vs_w501.txt"), sep='\t', index=True)
combo_counts_ISO1_vs_w501.to_csv(os.path.join(output_dir, "combo_counts_ISO1_vs_w501.txt"), sep='\t', index=True)
combo_counts_ISO1_vs_A4.to_csv(os.path.join(output_dir, "combo_counts_ISO1_vs_A4.txt"), sep='\t', index=True)


print("Script completed successfully.")
