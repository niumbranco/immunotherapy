# 🧬 16S Amplicon Sequencing Workflow (DADA2 + Snakemake)

[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![DADA2](https://img.shields.io/badge/DADA2-High%20Resolution%20Denoising-blue.svg)](https://benjjneb.github.io/dada2/index.html)

---

## Overview

This repository contains a workflow implemented in **Snakemake** to run **16S rRNA amplicon analyses**.  
It is designed for reproducible and automated microbiome workflows.

The pipeline performs:

- Read quality control with **fastp**
- Amplicon sequence variant (ASV) inference and taxonomic assignment with **DADA2**, using the **SILVA database [1]**
- Diversity and community composition visualizations

---

## Requirements

Before running the workflow, ensure that you have:

1. Raw 16S paired-end FASTQ files.
2. A table listing sample IDs and file names.
3. Snakemake installed and available in your system PATH.

Snakemake is the workflow manager that executes the pipeline and automatically creates all rule-specific environments (for `fastp`, `DADA2`, and R). It must be installed once before running the workflow. The recommended method is to create a dedicated Conda environment for it:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba -y
conda activate snakemake
snakemake --version
```
This ensures that Snakemake and its dependencies are correctly configured on your system.
You should see a version number (e.g., 7.32.4) when running the last command, confirming that the installation was successful.

## Repository Structure
```
16S_pipeline/
├── Snakefile                 # main Snakemake workflow
├── config.yml                # configuration file (paths, parameters)
├── envs/                     # conda environments per step
│   ├── fastp_env.yml
│   ├── dada2_env.yml
│   └── r_env.yml
├── database/                 # taxonomy database (e.g. SILVA)
│   ├── silva_nr99_v138.1_wSpecies_train_set.fa
│   └── silva_species_assignment_v138.1.fa
├── rules/                    # modular rule files
│   ├── fastp.rule
│   ├── dada2.rule
│   └── diversity_analysis.rule
├── scripts/                  # R scripts called by Snakemake
│   └── dada2_automated/
│       └── for_snakemake/
│           ├── prepare_tables_smk.R
│           ├── plot_taxonomy_composition_smk.R
│           ├── track_reads_smk.R
│           ├── assign_taxonomy_smk.R
│           └── (other DADA2 steps)
├── data/
│   ├── fastq/
│   │   ├── raw/              # input FASTQs
│   │   └── filtered/         # filtered reads after fastp
│   ├── fastp/                # QC reports (HTML/JSON)
│   └── reads_taxonomy/       # intermediate DADA2 outputs (seqtab, taxa, etc.)
├── results/                  # final outputs
│   ├── taxonomy_counts_P.tsv
│   ├── taxonomy_counts_G.tsv
│   └── visuals/
│       ├── microbial_composition_per_sample_P.pdf
│       ├── microbial_composition_overview_P.pdf
│       ├── microbial_composition_per_sample_G.pdf
│       └── microbial_composition_overview_G.pdf
└── README.md
```
---

### Conda environments
All environments to run the snakemake workflow can be built automatically using the provided `.yml` files within envs folder, or manually as described below.

#### 1. fastp environment

Used for quality filtering and integrity checks of FASTQ files.

```bash
conda create -n fastp_env -c bioconda fastp -y
conda activate fastp_env
fastp --version
```

#### 2. DADA2 environment

Used for sequence processing, ASV inference, and taxonomy assignment in R.
```bash
conda create -n r_env -c conda-forge r-base r-essentials -y
conda activate r_env
conda install -c bioconda -c conda-forge bioconductor-dada2 -y
```
This installs DADA2 and all dependencies directly into the R environment, ensuring compatibility with the R version and avoiding issues with install.packages() inside R.
To verify installation:
```bash
library(dada2)
packageVersion("dada2")
```

#### 3. Snakemake environment
Used to execute and manage the automated workflow.
```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake -y
conda activate snakemake
snakemake --version
```

### Results 
File	Description
taxonomy_counts_P.tsv	ASV counts per sample (Phylum level)
taxonomy_counts_G.tsv	ASV counts per sample (Genus level)
microbial_composition_per_sample_*.pdf	Stacked barplot per sample
microbial_composition_overview_*.pdf	Average composition across all samples

🧠 Notes
The workflow is modular — you can add or remove rule files easily (e.g., alpha_diversity.rule, differential_abundance.rule).
Each R script is executed through Snakemake, ensuring reproducible environments and file paths.
FASTQ files, intermediate .rds, and .snakemake/ folders should be excluded from git using .gitignore.
You can customize the plotting script to add sample grouping, metadata coloring, or horizontal bars.

👩‍🔬 Citations
[1] McLaren, M. R., & Callahan, B. J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4587955


