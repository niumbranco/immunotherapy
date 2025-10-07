import os
import pandas as pd
import csv

# Load the YAML into the config dictionary
configfile: "config/config.yml"

# Always relative to where the Snakefile lives
BASEDIR = workflow.basedir
ENV_DIR = os.path.join(BASEDIR, config["env_dir"])

# Definition of variables from config
input_dir = config['input_dir']
output_dir = config['output_dir']
env_dir = config['env_dir']
dada2_dir = config['dada2_dir']

sample_table = pd.read_csv(config["sample_table"], sep="\t", comment="#")

SAMPLES = sample_table["sample"].tolist()

# Function to generate input file paths based on sample ID
def get_input_files(wc):
    row = sample_table[sample_table["sample"] == wc.sample].iloc[0]
    return {"read_1": row["R1"], "read_2": row["R2"]}


# Generate list of expected outputs
expected_outputs = (
    expand(f"{output_dir}/{{sample}}_1.clean.fastq.gz", sample=SAMPLES) +
    expand(f"{output_dir}/{{sample}}_2.clean.fastq.gz", sample=SAMPLES) +
    expand(f"{dada2_dir}/{{sample}}_filt1.fastq.gz", sample=SAMPLES) +
    expand(f"{dada2_dir}/{{sample}}_filt2.fastq.gz", sample=SAMPLES) +
    [f"{dada2_dir}/errF.rds", f"{dada2_dir}/errR.rds"] +
    expand(f"{dada2_dir}/{{sample}}_dadaFs.rds", sample=SAMPLES) +
    expand(f"{dada2_dir}/{{sample}}_dadaRs.rds", sample=SAMPLES) +
    expand(f"{dada2_dir}/{{sample}}_mergers.rds", sample=SAMPLES) +
    [f"{dada2_dir}/seqtab.rds"] +
    [f"{dada2_dir}/taxa.rds"]

    #expand(f"{dada2_dir}/{{sample}}_table.rds", sample=SAMPLES) +
    #expand(f"{dada2_dir}/{{sample}}_seqs.rds", sample=SAMPLES)
)

print(f"fastp conda env: {env_dir}/fastp_env.yml")

#workdir:
#    output_dir

include: "rules/fastp.rule"
include: "rules/dada2.rule"
include: "rules/diversity_analysis.rule"

rule all:
     input:
        expected_outputs
