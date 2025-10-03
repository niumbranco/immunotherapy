#!/bin/bash

# Input metadata file
TSV="filereport_16S.tsv"

# Output log files
LOG_OK="md5_check_ok.txt"
LOG_FAIL="md5_check_fail.txt"

# Clear previous log files if they exist
> "$LOG_OK"
> "$LOG_FAIL"

echo "Checking .fastq.gz files using md5"

# Skip the header and process each line
tail -n +2 "$TSV" | while IFS=$'\t' read -r run_id sample exp study taxid organism md5s fastq_ftp _ _ _; do

    # Split URLs and MD5s (handle single or paired-end)
    IFS=';' read -r url1 url2 <<< "$fastq_ftp"
    IFS=';' read -r md5_1 md5_2 <<< "$md5s"

    file1=$(basename "$url1")
    file2=$(basename "$url2")

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
