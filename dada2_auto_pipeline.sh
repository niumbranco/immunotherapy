#!/bin/bash

# Define paths
root_dir="/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE"
read_trim_dir="${root_dir}/CODES/fastp_output"
results_dir="${root_dir}/CODES/dada2_output"
script_dir="${root_dir}/CODES/scripts/dada2_automated"
db_dir="${root_dir}/database"

# Create output directory if it doesn't exist
mkdir -p "$results_dir"

# Run filtering step
Rscript filt.R -i "$read_trim_dir" -o "$results_dir"

# Run error learning 
Rscript err.R -i "$results_dir" -o "$results_dir" 

# Run denoising step (ASV inference)
Rscript infer_ASV.R -i "$results_dir" -o "$results_dir"

while [ ! -f "$results_dir/dadaRs.rds" ]; do
  echo "Waiting for dadaRs.rds to be written..."
  sleep 2
done

# Construct the sequence table and remove chimeras 
Rscript seqtb.R -i "$read_trim_dir" -o "$results_dir"

# Taxonomic assignment with SILVA
Rscript "${script_dir}/silva.R" \
  -i "$results_dir" \
  -o "$results_dir" \
  -r "${db_dir}/silva_nr99_v138.1_wSpecies_train_set.fa.gz" \
  -s "${db_dir}/silva_species_assignment_v138.1.fa.gz"

