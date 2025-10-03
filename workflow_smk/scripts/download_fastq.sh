#!/bin/bash
while IFS= read -r url; do
	wget -c "$url"
done < fastq_links_16S.txt

