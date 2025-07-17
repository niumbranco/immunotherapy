import os
import pandas as pd

PWD = os.getcwd()

# Definition of environmental variables: paths for the source codes, among others
CONFIG = os.environ.get("config_file", "../config/config.yml")
configfile: CONFIG

# Definition of variables from config
tmp_dir = config['tmp_dir']
input_dir = config['input_dir']
output_dir = config['output_dir']
env_dir = config['env_dir']

sample_table = pd.read_csv(config["sample_table"], sep="\t", comment="#")

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
expected_outputs = []
for _, row in sample_table.iterrows():
    sample_id = row['sample']
    # Add trimmed fastq output files
    expected_outputs.append(f"{output_dir}/{sample_id}_R1.processed.fastq.gz")
    expected_outputs.append(f"{output_dir}/{sample_id}_R2.processed.fastq.gz")

print(f"fastp conda env: {env_dir}/fastp_env.yml")

workdir:
    output_dir

include: '../rules/run_trimming.smk'

rule all:
     input:
        expected_outputs
