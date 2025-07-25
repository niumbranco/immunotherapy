import os
import pandas as pd

PWD = os.getcwd()

# Definition of environmental variables: paths for the source codes, among others
CONFIG = os.environ.get("config_file", "../config/config.yml")
configfile: CONFIG

# Definition of variables from config
input_dir = config['input_dir']
output_dir = config['output_dir']
env_dir = config['env_dir']
dada2_dir = config['dada2_dir']

sample_table = pd.read_csv(config["sample_table"], sep="\t", comment="#")

SAMPLES = sample_table["sample"].tolist()

# Function to generate input file paths based on sample ID
def get_input_files(wildcards):
    sample_id = wildcards.sample
    try:
        row = sample_table[sample_table["sample"] == sample_id].iloc[0]
    except IndexError:
        raise ValueError(f"[ERROR] Sample '{sample_id}' not found in sample_table")
    try:
        read_1 = os.path.join(input_dir, row["R1"])
        read_2 = os.path.join(input_dir, row["R2"])
    except KeyError as e:
        raise KeyError(f"[ERROR] Column missing in sample_table: {e}")
    return {"read_1": read_1, "read_2": read_2}


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

workdir:
    output_dir

include: "rules/fastp.rule"
include: "rules/dada2.rule"

rule all:
     input:
        expected_outputs
