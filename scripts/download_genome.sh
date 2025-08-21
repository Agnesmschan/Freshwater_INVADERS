#!/bin/bash
#SBATCH --job-name=download_genome_fastas
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --array=0-9%5
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -euo pipefail

# ---- paths (adjust if needed) ----
ACC_LIST=$(readlink -f ../data/genome_accessions_to_download.txt)
DB_DIR=$(readlink -f ../kraken2_databases/Freshwater_invaders)
FA_DIR="$DB_DIR/fastas"
TMP_ROOT="$DB_DIR/tmp"
mkdir -p "$FA_DIR" "$TMP_ROOT"

# Load accessions and select the one for this task
mapfile -t ACCESSIONS < "$ACC_LIST"
i=${SLURM_ARRAY_TASK_ID}
if [ "$i" -ge "${#ACCESSIONS[@]}" ]; then
  echo "[ERROR] index $i >= ${#ACCESSIONS[@]}"; exit 1
fi
ACC="${ACCESSIONS[$i]}"
ACC="${ACC%%$'\r'}"   # strip CR if copied from Windows

echo "[INFO] task=$i accession=$ACC"

# Skip if we already have the fasta
OUT_FASTA="$FA_DIR/${ACC}.fna"
if [ -s "$OUT_FASTA" ]; then
  echo "[INFO] already exists: $OUT_FASTA (skipping)"
  exit 0
fi

# Per-task temp dir
TMP="$TMP_ROOT/temp_${SLURM_ARRAY_TASK_ID}"
rm -rf "$TMP"; mkdir -p "$TMP"

# Download + unzip
if ! datasets download genome accession "$ACC" --include genome --filename "$TMP/${ACC}.zip"; then
  echo "[WARN] download failed for $ACC"
  rm -rf "$TMP"
  exit 2
fi
unzip -q "$TMP/${ACC}.zip" -d "$TMP/${ACC}/"

# Concatenate all contigs into one fasta (handles .fna and .fna.gz)
touch "$OUT_FASTA"
find "$TMP/${ACC}/ncbi_dataset/data" -type f \( -name '*.fna' -o -name '*.fna.gz' \) -print0 \
| while IFS= read -r -d '' f; do
    case "$f" in
      *.gz) zcat "$f" >> "$OUT_FASTA" ;;
      *)    cat  "$f" >> "$OUT_FASTA" ;;
    esac
  done

if [ ! -s "$OUT_FASTA" ]; then
  echo "[WARN] no contigs found for $ACC"
  rm -f "$OUT_FASTA"
  rm -rf "$TMP"
  exit 3
fi

echo "[OK] wrote $OUT_FASTA ($(wc -c < "$OUT_FASTA") bytes)"
rm -rf "$TMP"
