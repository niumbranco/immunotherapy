
conda activate fastp_env

root_dir="/Users/smartinezarbas/Repositories/nium/roquette_nutriose"
data_dir="${root_dir}/data/raw_data"
read_report_dir="${root_dir}/data/fastp_reports"
read_trim_dir="${root_dir}/data/trimmed_reads"
results_dir="${root_dir}/data/results"
db_dir="/Users/smartinezarbas/databases"

mkdir -p ${read_report_dir} ${read_trim_dir} ${results_dir}

for sample in ${data_dir}/*_R1.fastq.gz
do
  base=$(basename ${sample} _R1.fastq.gz)

  fastp \
    -i ${data_dir}/${base}_R1.fastq.gz \
    -I ${data_dir}/${base}_R2.fastq.gz \
    -o ${read_trim_dir}/${base}_R1_trimmed.fastq.gz \
    -O ${read_trim_dir}/${base}_R2_trimmed.fastq.gz \
    -h ${read_report_dir}/${base}.html \
    -j ${read_report_dir}/${base}.json \
    --detect_adapter_for_pe \
    --length_required 150 \
    --cut_mean_quality 20 \
    --thread 4
done

conda deactivate

# conda activate r_env # with all the required packages 

# filter reads
Rscript dada2_filt.R -i ${read_trim_dir} \
                      -o ${results_dir} 

# learn errors
Rscript dada2_err.R -i ${read_trim_dir} \
                      -o ${results_dir} \
                      --silva_train_set ${db_dir}/silva_nr99_v138.1_train_set.fa

# get seqtable
Rscript dada2_derep.R -i ${read_trim_dir} \
                      -o ${results_dir} 

# get seqtb 
Rscript dada2_seqtb.R -i ${read_trim_dir} \
                      -o ${results_dir} 

# assign taxonomy
Rscript dada2_silva.R -i ${read_trim_dir} \
                      -o ${results_dir} \
                      --silva_train_set ${db_dir}/silva_nr99_v138.1_train_set.fa.gz \
                      --silva_species_set ${db_dir}/silva_species_assignment_v138.1.fa.gz
