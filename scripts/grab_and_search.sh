#!/usr/bin/env bash
#SBATCH --partition=medium
#SBATCH --output=sketch_and_gather.%A_%a.out
#SBATCH --error=sketch_and_gather.%A_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-12%5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

set -euo pipefail

# --- paths ---
FASTQ_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/sampled_eDNA_seq"
SIG_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/sampled_signatures"
OUT_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/sampled_smash_match"
TEMP_DIR="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/temporary_files"
DB1="/mnt/shared/projects/nhm/clark-student/MastersProjects/Summer-2025/Freshwater/data/Invaders_k31_database.sbt.zip"  # custom database

mkdir -p "$SIG_DIR" "$OUT_DIR" "$TEMP_DIR"

# --- params ---
K=31
SCALED=1000
THRESH_BP=1000

pad2() { printf "%02d" "$1"; }

# Collect FASTQs and sort *naturally* so filtered_10 > filtered_2
mapfile -t FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | sort -V)

if (( ${#FILES[@]} == 0 )); then
  echo "No FASTQ files found in $FASTQ_DIR"
  exit 1
fi

# Build stable per-barcode numbering so each barcodeNN gets 01,02,03...
declare -A bc_counts           # key=NN (digits), val=count so far
declare -A sample_idx_for_path # key=full path, val=per-barcode index (1-based)

for f in "${FILES[@]}"; do
  bn="$(basename "$f")"
  if [[ "$bn" =~ barcode([0-9]+) ]]; then
    bc="${BASH_REMATCH[1]}"
  else
    echo "WARNING: couldn't find 'barcodeNN' in: $bn; skipping numbering for this file."
    bc="0"
  fi
  n="${bc_counts[$bc]:-0}"
  n=$((n+1))
  bc_counts["$bc"]="$n"
  sample_idx_for_path["$f"]="$n"
done

# Bounds check for array index
if (( SLURM_ARRAY_TASK_ID >= ${#FILES[@]} )); then
  echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} >= number of files ${#FILES[@]}; exiting."
  exit 0
fi

FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
BN="$(basename "$FILE")"

# Extract barcode digits again for our file
if [[ "$BN" =~ barcode([0-9]+) ]]; then
  BC_NUM="${BASH_REMATCH[1]}"    # e.g., 57
  BC_LABEL="barcode${BC_NUM}"    # e.g., barcode57
else
  BC_NUM="0"
  BC_LABEL="barcode00"
fi

SNUM="${sample_idx_for_path["$FILE"]}"   # 1-based index within this barcode
SNUM_PAD="$(pad2 "$SNUM")"               # 01, 02, ...
BASENAME="${BC_LABEL}_${SNUM_PAD}"       # e.g., barcode57_01

OUT_SIG="${SIG_DIR}/${BASENAME}.sig"
OUTPUT_DATA="${TEMP_DIR}/${BASENAME}_data.csv"
FINAL_OUTPUT="${OUT_DIR}/${BASENAME}_results.csv"

echo "[$(date)] Processing: $BN -> ${BASENAME}"

# Create sourmash signature for this single FASTQ
# Use --name so the signature carries the same label
if [[ ! -s "$OUT_SIG" ]]; then
  sourmash sketch dna -p "scaled=${SCALED},k=${K},abund" --name "$BASENAME" -o "$OUT_SIG" "$FILE"
else
  echo "Signature exists, skipping sketch: $OUT_SIG"
fi

# Run gather against your custom DB
sourmash gather -k "$K" --threshold-bp "$THRESH_BP" "$OUT_SIG" "$DB1" -o "$OUTPUT_DATA"

# Prepend 'sample' column to the gather CSV
if [[ -s "$OUTPUT_DATA" ]]; then
  HEADER="$(head -n1 "$OUTPUT_DATA")"
  echo "sample,${HEADER}" > "$FINAL_OUTPUT"
  tail -n +2 "$OUTPUT_DATA" | awk -v sample="$BASENAME" -F, 'BEGIN{OFS=","} {print sample, $0}' >> "$FINAL_OUTPUT"
else
  echo "No rows from gather; writing header only."
  echo "sample" > "$FINAL_OUTPUT"
fi

rm -f "$OUTPUT_DATA"
echo "[$(date)] Done: $FINAL_OUTPUT"
