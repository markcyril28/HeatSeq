"""
matrix_builder.py

Description:
This script builds a gene expression matrix from multiple sample files.
- It takes a reference gene name listand several sample files (TSV, gene_name<TAB>count).
- For each gene in the reference, it finds the corresponding count in each sample file.
- If a gene is missing in a sample, '0' is inserted.
- Output is a tab-separated matrix: rows are genes, columns are samples.

Usage:
python matrix_builder.py gene_names.txt sample1.txt sample2.txt ...
"""

import sys

def main(gene_names_file, sample_files):
    # Read gene names
    with open(gene_names_file) as f:
        gene_names = [line.strip() for line in f]

    # Build a dict for each sample: {gene_name: count}
    sample_dicts = []
    for sample_file in sample_files:
        gene_to_count = {}
        with open(sample_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                gene = parts[0]
                count = parts[1]
                gene_to_count[gene] = count
        sample_dicts.append(gene_to_count)

    # Output matrix
    for gene in gene_names:
        row = [gene]
        for sample in sample_dicts:
            count = sample.get(gene, "0")
            row.append(count)
        print("\t".join(row))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python matrix_builder.py gene_names.txt sample1.txt sample2.txt ...")
        sys.exit(1)
    gene_names_file = sys.argv[1]
    sample_files = sys.argv[2:]
    main(gene_names_file, sample_files)
