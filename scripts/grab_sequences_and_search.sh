#!/bin/bash

#SBATCH --partition=medium
#SBATCH --output=grab_sequences_and_search.out
#SBATCH --error=grab_sequences_and_search.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-12%5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

FASTQ_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/sampled_eDNA_seq"
SIG_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/tmp_sampled_signatures"
OUT_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/sampled_smash_match"
TEMP_DIR="../data/temporary_files"
DB1="../data/Invaders_k31_database.sbt.zip" #Custom database

mkdir -p "$SIG_DIR" "$OUT_DIR" "$TEMP_DIR"

BARCODE=

# Create sourmash signature
sourmash sketch dna -p scaled=1000,k=31,abund "$FASTQ_DIR"/*.fastq -o "$SIG_DIR/${BARCODE}.sig"


TEMP_DIR="../data/temporary_files"
DB1="../data/Invaders_k31_database.sbt.zip" #Custom database

# Gather all .sig files into an array
mapfile -t FILES < <(find "$SIG_DIR" -maxdepth 1 -name '*.sig' | sort)

# Check files...
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#FILES[@]}" ]; then
  echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} exceeds number of input files (${#FILES[@]}; aborting..)"
  exit 1
fi

FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$FILE" .sig)

# Output file paths
OUTPUT_DATA="${TEMP_DIR}/${BASENAME}_data.csv"
FINAL_OUTPUT="${OUT_DIR}/${BASENAME}_results.csv"

echo "[$(date)] Starting sourmash gather for: $BASENAME"

# Run search with sourmash gather
sourmash gather -k 31 --threshold-bp 1000 "$FILE" "$DB1" -o "$OUTPUT_DATA"

# Read the header from OUTPUT_DATA, and prepend sample column
HEADER=$(head -n1 "$OUTPUT_DATA")


echo "sample,${HEADER}" > "$FINAL_OUTPUT"

# For each data row, prepend the sample name ($BASENAME) and skip header
tail -n +2 "$OUTPUT_DATA" | \
  awk -v sample="$BASENAME" -F, 'BEGIN{OFS=","} {print sample, $0}' \
  >> "$FINAL_OUTPUT"

# 3. Clean up
rm -f "$OUTPUT_DATA"

echo "[$(date)] Finished processing $BASENAME"
