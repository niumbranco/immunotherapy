#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Utility for Snakemake workflows. Database download
# -----------------------------------------------------------------------------
# This script downloads files from provided URLs,
# decompresses them if necessary, and moves them into the local "database/"
# directory to run the workflow. File names are automatically inferred from the URLs.
#
# Usage:
#   bash scripts/download_database.sh <URL1> [URL2] [URL3] ...
#
# Example:
#   bash scripts/download_database.sh \
#     https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz \
#     https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
# -----------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <URL1> [URL2] ..."
    exit 1
fi

# Define destination directory
DB_DIR="database"
mkdir -p "${DB_DIR}"

# Check if wget or curl is available
if command -v wget >/dev/null 2>&1; then
    DOWNLOADER="wget -q -O"
elif command -v curl >/dev/null 2>&1; then
    DOWNLOADER="curl -s -L -o"
else
    echo "Error: neither wget nor curl found. Please install one and rerun." >&2
    exit 1
fi

echo "----------------------------------------------------------"
echo "Downloading database files..."
echo "Destination: ${DB_DIR}"
echo "----------------------------------------------------------"

for URL in "$@"; do
    FILENAME=$(basename "${URL}")
    DEST="${DB_DIR}/${FILENAME}"

    echo "→ Downloading ${FILENAME}"
    ${DOWNLOADER} "${DEST}" "${URL}"

    # Decompress if gzipped
    if [[ "${DEST}" == *.gz ]]; then
        echo "  Decompressing ${FILENAME}"
        gunzip -f "${DEST}"
    fi
done

echo "----------------------------------------------------------"
echo "All downloads completed. Contents of ${DB_DIR}:"
ls -lh "${DB_DIR}"
echo "----------------------------------------------------------"
