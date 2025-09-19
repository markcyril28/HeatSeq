#!/bin/bash 

# Gunzipped Files Copier. 
scp -r \
    admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/*.gz \
    /mnt/f/MCRM_pipeline_transient_folder

scp -r \
    admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/*.gz \
    /mnt/f/MCRM_pipeline_transient_folder
