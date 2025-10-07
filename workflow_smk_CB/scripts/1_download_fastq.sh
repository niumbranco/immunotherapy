#!/bin/bash

## 1. Fetch data
# Create directory for output
OUTDIR="../data/fastq"
LOGFILE="$OUTDIR/download_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$OUTDIR"

echo "Download started at $(date)" > "$LOGFILE"
echo "Output directory: $OUTDIR" >> "$LOGFILE"
echo "---------------------------------" >> "$LOGFILE"

# Read URLs and download
while IFS= read -r url; do
    echo "Processing $url ..." | tee -a "$LOGFILE"
    if wget -c -P "$OUTDIR" "$url" >> "$LOGFILE" 2>&1; then
        echo "SUCCESS: $(basename "$url") downloaded" | tee -a "$LOGFILE"
    else
        echo "FAILED: $(basename "$url")" | tee -a "$LOGFILE"
    fi
    echo "---------------------------------" >> "$LOGFILE"
done < ../data/fastq_links_16S.txt

echo "Download finished at $(date)" >> "$LOGFILE"

## 2. Check integrity
# Input metadata file
TSV="../data/filereport_16S.tsv" # this file came, in this case, from the metadata associated with the project

# The right location of the fastq files
FASTQ_DIR="../data/fastq"

# Output log files
LOG_OK="md5_check_ok.txt"
LOG_FAIL="md5_check_fail.txt"

# Clear previous log files if they exist
: > "$LOG_OK"
: > "$LOG_FAIL"

echo "Checking .fastq.gz files using md5"

# Skip the header and process each line
tail -n +2 "$TSV" | while IFS=$'\t' read -r run_id sample exp study taxid organism md5s fastq_ftp _ _ _; do

    # Split URLs and MD5s (handle single or paired-end)
    IFS=';' read -r url1 url2 <<< "$fastq_ftp"
    IFS=';' read -r md5_1 md5_2 <<< "$md5s"

    file1="$FASTQ_DIR/$(basename "$url1")"
    file2="$FASTQ_DIR/$(basename "$url2")"

    # Check file1
    if [[ -f "$file1" ]]; then
        md5_local=$(md5 -q "$file1")
        if [[ "$md5_local" == "$md5_1" ]]; then
            echo "$file1 OK" >> "$LOG_OK"
        else
            echo "$file1 FAILED (expected $md5_1, got $md5_local)" >> "$LOG_FAIL"
        fi
    else
        echo "$file1 MISSING" >> "$LOG_FAIL"
    fi

    # Check file2 (if exists)
    if [[ -n "$file2" ]]; then
        if [[ -f "$file2" ]]; then
            md5_local=$(md5 -q "$file2")
            if [[ "$md5_local" == "$md5_2" ]]; then
                echo "$file2 OK" >> "$LOG_OK"
            else
                echo "$file2 FAILED (expected $md5_2, got $md5_local)" >> "$LOG_FAIL"
            fi
        else
            echo "$file2 MISSING" >> "$LOG_FAIL"
        fi
    fi

done

echo "Check complete. See $LOG_OK and $LOG_FAIL for results."
