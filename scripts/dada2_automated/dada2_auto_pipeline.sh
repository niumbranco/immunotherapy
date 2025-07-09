#!/bin/bash

# Define paths
root_dir="/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE"
read_trim_dir="${root_dir}/CODES/fastp_output"
results_dir="${root_dir}/CODES/dada2_output"
script_dir="${root_dir}/CODES/scripts/dada2_automated"

# Create output directory if it doesn't exist
mkdir -p "$results_dir"

# Run filtering step
Rscript filt.R -i "$read_trim_dir" -o "$results_dir"
