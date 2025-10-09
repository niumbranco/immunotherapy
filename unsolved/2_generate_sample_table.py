import pandas as pd
import os

# Entrées
metadata_file = "/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/metadata/filereport_16S.tsv"
reads_dir = "/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/PRJEB61942/raw_reads"
output_file = "/Users/Emma/Documents/EIDD/2A/STAGE_2A/DOCUMENTS_DE_STAGE/immunotherapy/PRJEB61942/samples_ID.tsv"

# Lecture des métadonnées
df = pd.read_csv(metadata_file, sep="\t")

# Construction de la table avec vérification de l'existence des fichiers
rows = []
for run in df["run_accession"]:
    r1 = os.path.join(reads_dir, f"{run}_1.fastq.gz")
    r2 = os.path.join(reads_dir, f"{run}_2.fastq.gz")
    if os.path.exists(r1) and os.path.exists(r2):
        rows.append({"sample": run, "R1": r1, "R2": r2})
    else:
        print(f"[SKIPPED] Fichiers manquants pour {run}")

# Export
sample_df = pd.DataFrame(rows)
sample_df.to_csv(output_file, sep="\t", index=False)
print(f"Table enregistrée dans : {output_file}")