#!/bin/bash

# Define paths for test
root_dir="/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE"
read_trim_dir="${root_dir}/CODES/fastp_output_test"     # dossier avec 4-5 fastq
results_dir="${root_dir}/CODES/dada2_output_test"
script_dir="${root_dir}/CODES/scripts/dada2_automated"

mkdir -p "$results_dir"

# Run only filtering step for now
Rscript "$script_dir/dada2_filt.R" \
  -i "$read_trim_dir" \
  -o "$results_dir"
