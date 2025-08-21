#!/bin/bash

#SBATCH --partition=medium
#SBATCH --output=preget_samples_p.out
#SBATCH --error=preget_samples_p.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1953%20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

set -euo pipefail

# Set up directory

ACCESSION_LIST="../data/Freshwater_accessions.txt"
mkdir -p ../data/test_sra_signatures
mkdir -p ./temp_files2 logs

# Set error parameters to terminate job when list not found

if [ ! -f "$ACCESSION_LIST" ]; then
   echo "ERROR: Accession list $ACCESSION_LIST not found."
   exit 1

else
   echo "Found $ACCESSION_LIST, now processing..."
fi

# Load accessions and assign the one for this task

mapfile -t ACCESSIONS < "$ACCESSION_LIST"
ACCESSION=${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

# Store in Temp directory
TEMP_DIR="./temp_files2/temp_$SLURM_ARRAY_TASK_ID"
mkdir -p "$TEMP_DIR"
    echo "[$(date)] Processing $ACCESSION"

# Start downloading SRA

if ! prefetch --max-size 500G "$ACCESSION" -o "$TEMP_DIR/${ACCESSION}.sra"; then
  echo "$ACCESSION download failed" >> logs/failed_downloads.txt
  rm -rf "$TEMP_DIR"
  exit 1
fi

# Convert to FASTQ
if ! fasterq-dump "$TEMP_DIR/${ACCESSION}.sra" --concatenate-reads --skip-technical --outdir "$TEMP_DIR"; then
  echo "$ACCESSION FASTQ conversion failed" >> logs/failed_fastq.txt
  rm -rf "$TEMP_DIR"
  exit 1
fi

# Create sourmash signature
sourmash sketch dna -p scaled=1000,k=31,abund "$TEMP_DIR"/*.fastq -o ../data/test_sra_signatures/"${ACCESSION}.sig"

echo "$ACCESSION processed successfully" >> logs/success.log
rm -rf "$TEMP_DIR"

