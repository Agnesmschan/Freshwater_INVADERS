#!/bin/bash

#SBATCH --partition=medium
#SBATCH --job-name=check_database
#SBATCH --output=check_database.out
#SBATCH --error=check_database.err
#SBATCH --mem-per-cpu=130G
#SBATCH --cpus-per-task=1

set -euo pipefail

DB="/mnt/shared/projects/nhm/clark-student/AirSeq/Reference_Databases/core_nt_database/core_nt_20250609"
TARGETS="../kraken2_databases/Freshwater_invaders/metadata/species_taxids.tsv"

# -- outputs (write to your space, not into DB) --
OUTDIR="../data/check_db_results"
mkdir -p "$OUTDIR"

# sanity checks
[[ -d "$DB" ]] || { echo "[ERROR] DB not found: $DB" >&2; exit 1; }
[[ -f "$DB/names.dmp" ]] || { echo "[ERROR] names.dmp missing in $DB" >&2; exit 1; }
[[ -f "$TARGETS" ]] || { echo "[ERROR] targets file not found: $TARGETS" >&2; exit 1; }

# A) Taxonomy presence check (use scientific names only from names.dmp; it's pipe-delimited)
awk -F '\t\\|\\t' '$4 ~ /scientific name/ {gsub(/^ +| +$/,"",$1); gsub(/^ +| +$/,"",$2); print $1"\t"$2 }' \
  "$DB/names.dmp" > "$OUTDIR/names_scientific.tsv"
# file: taxid \t scientific_name

echo "[INFO] Checking taxonomy presence for targets…"
awk -F'\t' 'NR>1 {print $3"\t"$2"\t"$1}' "$TARGETS" | while IFS=$'\t' read -r tid species acc; do
  hit=$(awk -F'\t' -v t="$tid" '$1==t {print $2; exit}' "$OUTDIR/names_scientific.tsv" || true)
  if [[ -n "$hit" ]]; then
    printf "TAXONOMY_OK\t%s\t%s\t%s\n" "$tid" "$species" "$acc"
  else
    printf "TAXONOMY_MISSING\t%s\t%s\t%s\n" "$tid" "$species" "$acc"
  fi
done > "$OUTDIR/taxonomy_check.tsv"

# B) Database content presence (does the DB actually have k-mers at/under that taxid?)
echo "[INFO] Running kraken2-inspect (this only reads the DB)…"
kraken2-inspect --db "$DB" --memory-mapping --report-zero-counts > "$OUTDIR/inspect_full.tsv"

# kraken2-inspect output columns (tab-separated):
# percent   clade_reads   taxon_reads   rank   taxid   name(with indentation)
# We will report PRESENT if clade_reads > 0 for the taxid.

# Build a compact summary for your targets
awk -F'\t' 'NR>1 {print $3"\t"$2"\t"$1}' "$TARGETS" > "$OUTDIR/targets.list"
# targets.list: taxid \t species \t accession

awk -F'\t' 'NR==FNR{meta[$1]=$2"\t"$3; next}
            {
              tid=$5; clade=$2; name=$6; gsub(/^ +/,"",name)
              if (tid in meta) {
                present = (clade+0)>0 ? "PRESENT" : "ZERO";
                print tid "\t" meta[tid] "\t" present "\t" clade "\t" name
              }
            }' "$OUTDIR/targets.list" "$OUTDIR/inspect_full.tsv" \
  | awk 'BEGIN{OFS="\t"; print "taxid","species","accession","status","clade_reads","db_name"}1' \
  > "$OUTDIR/targets_in_db.tsv"

echo "[INFO] Wrote:"
echo " - $OUTDIR/taxonomy_check.tsv"
echo " - $OUTDIR/targets_in_db.tsv"
echo " - $OUTDIR/inspect_full.tsv (full inspect report)"

# Quick preview
echo "----- taxonomy_check (head) -----"
head "$OUTDIR/taxonomy_check.tsv" || true
echo "----- targets_in_db (head) -----"
head "$OUTDIR/targets_in_db.tsv" || true
