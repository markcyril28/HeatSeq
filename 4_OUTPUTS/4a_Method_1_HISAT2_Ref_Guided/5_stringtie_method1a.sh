#!/bin/bash

# Reference Guided: 
# Transcript Assembly and Quantification (StringTie)
stringtie sample.bam -G reference.gtf -o sample.gtf -A sample.gene_abund.tab

# Generate Expression Matrix (prepDE.py)
python prepDE.py -i sample_list.txt

: << 'EXAMPLE'
sample1 sample1.gtf
sample2 sample2.gtf
...
EXAMPLE
