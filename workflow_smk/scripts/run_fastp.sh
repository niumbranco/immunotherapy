#!/bin/bash

#create output and reprot folders
mkdir -p fastp_output 
mkdir -p fastp_reports 

#loop over all _1.fastq.gz files 
for file1 in *_1.fastq.gz; do 
	sample=${file1%%_1.fastq.gz} #get the base name 
	file2="${sample}_2.fastq.gz" #define corresponding file2
	
	#define output files
	out1="fastp_output/${sample}_1.clean.fastq.gz"
	out2="fastp_output/${sample}_2.clean.fastq.gz"

	#define report files
	html_report="fastp_reports/${sample}_report.html"
	json_report="fastp_reports/${sample}_report.json"

	#run fastp
	fastp -i "$file1" -I "$file2" -o "$out1" -O "$out2" \
		-h "$html_report" -j "$json_report" \
		--thread 4
	echo "Done with $sample"
done
