import csv

# Define a dictionary to store ortholog pairs
ortholog_dict = {}

# Read the ortholog pairs from the second file
with open('mel-sim_orthologs.txt', 'r') as ortholog_file:
    for line in ortholog_file:
        line = line.strip().split('\t')
        ortholog_dict[line[1]] = line[0]

# Create a list to store rows with orthologs
ortholog_rows = []

# Counters to track replaced genes and missing gene-ortholog pairs
replaced_genes_count = 0
missing_pairs_count = 0

# Read the first file and replace gene names with orthologs
with open('all_genes_w501.csv', 'r') as gene_file:
    csv_reader = csv.reader(gene_file)
    header = next(csv_reader)  # Assumes file has a header
    ortholog_rows.append(header)  # Retain header

    for row in csv_reader:
        gene_name = row[1]
        if gene_name in ortholog_dict:
            row[1] = ortholog_dict[gene_name]
            replaced_genes_count += 1
        else:
            missing_pairs_count += 1
        ortholog_rows.append(row)

# Write the ortholog rows to a new CSV file
with open('all_genes_orthologs.csv', 'w', newline='') as ortholog_csv:
    csv_writer = csv.writer(ortholog_csv)
    csv_writer.writerows(ortholog_rows)

# Print the count of replaced genes and missing gene-ortholog pairs
print("Number of genes replaced with orthologs:", replaced_genes_count)
print("Number of gene-ortholog pairs not found in all_genes.csv:", missing_pairs_count)
