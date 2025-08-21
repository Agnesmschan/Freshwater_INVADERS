#!/bin/bash
#SBATCH --job-name=kraken_class_array
#SBATCH --partition=medium
#SBATCH --cpus-per-task=6
#SBATCH --mem=256G
#SBATCH --output=kraken2_core.out
#SBATCH --error=kraken2_core.err
#SBATCH --array=0-1953%1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mansum.chan@nhm.ac.uk

set -euo pipefail

# Set directory and paths
ACCESSION_LIST="../data/Freshwater_accessions.txt"
TARGETS_TSV="../kraken2_databases/Freshwater_invaders/metadata/species_taxids.tsv"
DB="/mnt/shared/projects/nhm/clark-student/AirSeq/Reference_Databases/core_nt_database/core_nt_20250609"

# Outputs
WORKROOT="../data"
FQ_DIR="$WORKROOT/sra_fastq_core"
TMP_BASE="$WORKROOT/tmp"
REPORT_DIR="$WORKROOT/kraken_reports_core"
KRAKEN_DIR="$WORKROOT/kraken_raw_core"
TARGETS_DIR="$WORKROOT/kraken_targets_core"
LOG_DIR="logs"

mkdir -p "$FQ_DIR" "$TMP_BASE" "$REPORT_DIR" "$KRAKEN_DIR" "$TARGETS_DIR" "$LOG_DIR"

# Validate inputs files
if [[ ! -f "$ACCESSION_LIST" ]]; then
  echo "[ERROR] Accession list not found: $ACCESSION_LIST" >&2
  exit 1
fi
if [[ ! -f "$TARGETS_TSV" ]]; then
  echo "[ERROR] targets.tsv not found: $TARGETS_TSV" >&2
  exit 1
fi

# Build targets.taxids once (safe to redo in each task; it’s cheap)
TARGETS_TAXIDS="$WORKROOT/targets.taxids"
awk -F'\t' 'NR>1 && $3!="" {print $3}' "$TARGETS_TSV" > "$TARGETS_TAXIDS"

# Get accession
mapfile -t ACCESSIONS < "$ACCESSION_LIST"
if (( SLURM_ARRAY_TASK_ID >= ${#ACCESSIONS[@]} )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID out of range (${SLURM_ARRAY_TASK_ID} >= ${#ACCESSIONS[@]})" >&2
  exit 2
fi
ACC="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"

# Make tem dir per-task
TMP="$TMP_BASE/${ACC}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMP"

echo "[$(date)] START $ACC  (task ${SLURM_ARRAY_TASK_ID})"

# 1) Download SRA & convert to FASTQ
echo "[INFO] prefetch $ACC"
prefetch --max-size 500G -O "$TMP" "$ACC"

SRA_FILE="$TMP/$ACC/$ACC.sra"
if [[ ! -s "$SRA_FILE" ]]; then
  # some sratoolkit versions place .sra directly under TMP
  SRA_FILE="$TMP/$ACC.sra"
fi
if [[ ! -s "$SRA_FILE" ]]; then
  echo "[ERROR] cannot find downloaded SRA for $ACC" >&2
  exit 3
fi

echo "[INFO] fasterq-dump $ACC"
fasterq-dump --split-files --skip-technical \
             -e "$SLURM_CPUS_PER_TASK" \
             -O "$TMP" "$SRA_FILE"

# compress FASTQs
if command -v pigz >/dev/null 2>&1; then
  pigz -p "$SLURM_CPUS_PER_TASK" "$TMP/${ACC}"_*.fastq
else
  find "$TMP" -maxdepth 1 -name "${ACC}_*.fastq" -exec gzip -1 {} +
fi

# move to final fastq dir
for f in "$TMP/${ACC}"_*.fastq.gz; do
  mv -f "$f" "$FQ_DIR/"
done

# identify R1/R2 (paired) or single-end
R1="$FQ_DIR/${ACC}_1.fastq.gz"
R2="$FQ_DIR/${ACC}_2.fastq.gz"
SE="$FQ_DIR/${ACC}.fastq.gz"

# 2) Kraken2 classification
# -------------------------------
REPORT="$REPORT_DIR/${ACC}.kreport"
KRAW="$KRAKEN_DIR/${ACC}.kraken"

if [[ -s "$R1" && -s "$R2" ]]; then
  echo "[INFO] kraken2 paired $ACC"
  kraken2 --db "$DB" \
          --memory-mapping \
          --threads "$SLURM_CPUS_PER_TASK" \
          --paired --gzip-compressed \
          --use-names \
          --confidence 0.1 \
          --report "$REPORT" \
          --output "$KRAW" \
          "$R1" "$R2"
elif [[ -s "$SE" ]]; then
  echo "[INFO] kraken2 single-end $ACC"
  kraken2 --db "$DB" \
          --memory-mapping \
          --threads "$SLURM_CPUS_PER_TASK" \
          --gzip-compressed \
          --use-names \
          --confidence 0.1 \
          --report "$REPORT" \
          --output "$KRAW" \
          "$SE"
else
  echo "[ERROR] no FASTQ found for $ACC in $FQ_DIR" >&2
  exit 4
fi

# 3) Filter to your target species & tidy CSV
TREP="$REPORT_DIR/${ACC}.targets.kreport"
CSV="$TARGETS_DIR/${ACC}.targets.csv"

# keep only lines whose taxid is in your targets file (column 5 in kreport)
awk -F'\t' 'NR==FNR{t[$1]; next} ($5 in t)' "$TARGETS_TAXIDS" "$REPORT" > "$TREP"

# If kreport is empty, delete and skip CSV creation
if [[ ! -s "$TREP" ]]; then
  echo "[INFO] No target species found for $ACC — removing empty $TREP"
  rm -f "$TREP"
else
  # make a clean CSV: sample, taxid, name, percent, clade_reads, direct_reads
  awk -F'\t' -v s="$ACC" 'BEGIN{OFS=","; print "sample,taxid,name,percent,clade_reads,direct_reads"}
                          {print s,$5,$6,$1,$2,$3}' "$TREP" > "$CSV"

  # If CSV is empty (only header), delete it
  if [[ $(wc -l < "$CSV") -le 1 ]]; then
    echo "[INFO] $CSV contains only header — removing"
    rm -f "$CSV"
  fi
fi

echo "[$(date)] DONE $ACC"

# 4) Cleanup temp
rm -rf "$TMP"

