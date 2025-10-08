🧬 16S Amplicon Sequencing Workflow (DADA2 + Snakemake)

This repository contains a **Snakemake-based 16S rRNA amplicon analysis pipeline**, designed for reproducible and automated microbiome workflows.  
It performs read quality control with **fastp**, denoising with **DADA2**, taxonomic assignment using the **SILVA database [1]**, and produces diversity visualizations.

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![DADA2](https://img.shields.io/badge/DADA2-%20v1.28%20|%20R%20package-blue.svg)](https://benjjneb.github.io/dada2/index.html)
---

## ⚙️ Pipeline Overview

| Step | Description | Tool |
|------|--------------|------|
| 1 | Quality control and read trimming | `fastp` |
| 2 | Denoising, merging, chimera removal and taxonomic assignment | `DADA2` |
| 3 | Diversity and community composition plots | `R` |

---

## 📁 Repository Structure
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

## 🧩 Installation

### 1️⃣ Clone the repository
```bash
git clone https://github.com/<yourusername>/16S_pipeline.git
cd 16S_pipeline
2️⃣ Create Conda environments
conda env create -f envs/fastp_env.yml
conda env create -f envs/dada2_env.yml
conda env create -f envs/r_env.yml
3️⃣ Set up database
Place the SILVA reference files in the database/ folder:
silva_nr99_v138.1_wSpecies_train_set.fa
silva_species_assignment_v138.1.fa
4️⃣ Configure the workflow
Edit config.yml to specify your paths and parameters.
Example:
fastp:
  quality_threshold: 20
  length_threshold: 100
  threads: 4

dada2:
  trunc_len_f: 250
  trunc_len_r: 200
  threads: 4

paths:
  raw_dir: "data/fastq/raw"
  output_dir: "results"
🚀 Running the Workflow
Run the full pipeline
snakemake --cores 4 --use-conda
Run only the plotting step
snakemake --cores 2 --use-conda results/visuals/microbial_composition_overview_G.pdf
All results will be saved in the results/ directory.

📊 Outputs
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

👩‍🔬 Citation
Callahan BJ et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods 13(7):581–583.
[1] McLaren, M. R., & Callahan, B. J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4587955
Köster J & Rahmann S (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics 28(19):2520–2522.


📄 License
This workflow is distributed under the MIT License.
You are free to reuse, modify, and redistribute it.
✨ Acknowledgements
Developed by Catarina Banco
Designed for automated 16S microbiome analyses at NIUM / Friesland Campina.