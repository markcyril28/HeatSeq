#!/bin/bash
# Method 3 Trinity - Conda Environment Setup
# Creates conda environment for Trinity de novo assembly

conda create -n trinity_env trinity hisat2 samtools stringtie salmon rsem r-base r-deseq2 bioconductor-tximport -y
