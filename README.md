# Immunotherapy Microbiome Analysis (PRJEB61942) 

This repository contains all scripts and ressources related to the analysis of the PRJEB61942 dataset, which investigates the link between the microbiome and response to immunotherapy in melanoma patients.


## Workflow: 

### 1. Data selection
From the ENA metadata, I selecting only the **16S rRNA** samples based on the "tax_id=9606" and confirmed the sequencing method 

### 2. Downloading the FASTQ files 
Using a script ('downlaod_fastq.sh'), I downloaded all paired-end FASTQ files related to 16S sequencing from the ENA FTP links

### 3. Files integrity check 
To make sure files were correctly downloaded: 
- I used the **md5 checksums** provided in the metadata
- I wrote a bash script ("verify_md5.sh") to compare them with local files
- A log of passed and failed checks was created
  
### 4. Quality control with fastp 
I ran fastp on all paired-end FASTQ files unsing default settings 
It cleaned low-quality reads and generated reports. 


## File structure: 

### scripts:
* download_fastq.sh    #download FATSQ files from ENA
* verify_md5.sh        #check the md5 integrity of downloaded files
* run_fastp.sh         #quality filtering with fastp

### metadata: 
* filereport_16S.txt   #cleaned metadata with 16S samples
* fastq_links_16S.txt  #list of FASTQ URLs to download

### logs: 
* md5_check_ok.txt     #files that passed the md5 check
* md5_check_fail.txt   #files that failed or were missing
* log_download.txt     #output log of the download script

### fastp      
* fastp_output        #cleaned FASTQ files
* fastp_reports       #HTML/JSON quality reports 
