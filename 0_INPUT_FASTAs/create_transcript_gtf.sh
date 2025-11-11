#!/bin/bash
# Create transcript-based GTF from FASTA for reference-guided quantification

FASTA="All_Smel_Genes_Full_Name.fasta"
OUTPUT="All_Smel_Genes_transcript_based.gtf"

awk '
BEGIN {OFS="\t"}
/^>/ {
    if (seq_name != "") {
        len = length(sequence)
        gene_id = gensub(/\.01$/, "", "g", seq_name)
        print seq_name, "custom", "transcript", 1, len, ".", "+", ".", "transcript_id \"" seq_name "\"; gene_id \"" gene_id "\";"
        print seq_name, "custom", "exon", 1, len, ".", "+", ".", "transcript_id \"" seq_name "\"; gene_id \"" gene_id "\";"
    }
    seq_name = substr($1, 2)
    sequence = ""
    next
}
{
    sequence = sequence $0
}
END {
    if (seq_name != "") {
        len = length(sequence)
        gene_id = gensub(/\.01$/, "", "g", seq_name)
        print seq_name, "custom", "transcript", 1, len, ".", "+", ".", "transcript_id \"" seq_name "\"; gene_id \"" gene_id "\";"
        print seq_name, "custom", "exon", 1, len, ".", "+", ".", "transcript_id \"" seq_name "\"; gene_id \"" gene_id "\";"
    }
}
' "$FASTA" > "$OUTPUT"

echo "Created: $OUTPUT"
wc -l "$OUTPUT"
