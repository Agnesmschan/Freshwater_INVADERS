#!/bin/bash

#SBATCH --job-name=final_search_run
#SBATCH --output=../data/temporary_files/slurm-%a.out
#SBATCH --error=../data/temporary_files/slurm-%a.err
#SBATCH --array=0-1943%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

set -euo pipefail

# Define input and output directories
SIG_DIR="../data/Freshwater_sra_signatures"
OUT_DIR="../data/outputs_final"
TEMP_DIR="../data/temporary_files"
DB1="../data/Invaders_k31_database.sbt.zip" #Custom database

mkdir -p "$OUT_DIR" "$TEMP_DIR"

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
