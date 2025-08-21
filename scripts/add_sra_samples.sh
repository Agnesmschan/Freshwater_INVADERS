#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name=add_sra_sketch
#SBATCH --output=add_sra_sketch.out
#SBATCH --error=add_ara_sketch.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-4764%20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

set -euo pipefail

# ---- config & paths ----
ACCESSION_LIST="../data/Freshwater_accessions_large.txt"
SIG_DIR="../data/Freshwater_sra_signatures"
WORK_TMP="./temp_files"
LOG_DIR="./logs"

mkdir -p "$SIG_DIR" "$WORK_TMP" "$LOG_DIR"

# ---- sanity ----
if [[ ! -f "$ACCESSION_LIST" ]]; then
  echo "[ERROR] Accession list not found: $ACCESSION_LIST" >&2
  exit 1
fi

# ---- pick accession for this task ----
mapfile -t ACCESSIONS < "$ACCESSION_LIST"
if (( SLURM_ARRAY_TASK_ID >= ${#ACCESSIONS[@]} )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} out of range (${#ACCESSIONS[@]} accessions)." >&2
  exit 2
fi
ACC="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"

# skip if signature already exists
OUT_SIG="$SIG_DIR/${ACC}.sig"
if [[ -s "$OUT_SIG" ]]; then
  echo "[$(date)] SKIP: $ACC (signature exists: $OUT_SIG)"
  exit 0
fi

# temp workspace
TMP="$WORK_TMP/temp_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

echo "[$(date)] START: $ACC (task ${SLURM_ARRAY_TASK_ID})"

# simple retry helper
retry() {
  local n=0 max=2 delay=20
  until "$@"; do
    n=$((n+1))
    if (( n > max )); then return 1; fi
    echo "[WARN] retry $n/$max: $*" >&2
    sleep "$delay"
  done
}

# 1) download SRA
# use -O <dir> and then resolve where the .sra landed
if ! retry prefetch --max-size 500G -O "$TMP" "$ACC"; then
  echo "$ACC download failed" >> "$LOG_DIR/failed_downloads.txt"
  exit 3
fi

SRA="$TMP/$ACC/$ACC.sra"
if [[ ! -s "$SRA" ]]; then
  SRA="$TMP/$ACC.sra"
fi
if [[ ! -s "$SRA" ]]; then
  # last resort: find it
  SRA=$(find "$TMP" -maxdepth 2 -type f -name "${ACC}.sra" -print -quit || true)
fi
if [[ -z "${SRA:-}" || ! -s "$SRA" ]]; then
  echo "[ERROR] cannot find downloaded SRA for $ACC under $TMP" >&2
  echo "$ACC SRA missing after prefetch" >> "$LOG_DIR/failed_downloads.txt"
  exit 3
fi

# 2) convert to FASTQ (concatenated reads; good for sourmash)
if ! retry fasterq-dump "$SRA" --concatenate-reads --skip-technical --outdir "$TMP"; then
  echo "$ACC FASTQ conversion failed" >> "$LOG_DIR/failed_fastq.txt"
  exit 4
fi

# 3) sketch with sourmash
shopt -s nullglob
FASTQS=("$TMP"/*.fastq)
if (( ${#FASTQS[@]} == 0 )); then
  echo "[ERROR] no FASTQ created for $ACC" >&2
  exit 5
fi

sourmash sketch dna -p scaled=1000,k=31,abund "${FASTQS[@]}" -o "$OUT_SIG"

rm -rf "$TMP"

echo "$ACC processed successfully" >> "$LOG_DIR/success.log"
echo "[$(date)] DONE: $ACC"
