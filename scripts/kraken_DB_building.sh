#!/bin/bash
#SBATCH --job-name=kraken_DB_building
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -euo pipefail

# ---- directory paths ----
ACC_LIST=$(readlink -f ../data/genome_accessions_to_download.txt)
DB_DIR=$(readlink -f ../kraken2_databases/Freshwater_invaders)
FA_DIR="$DB_DIR/fastas"
META="$DB_DIR/metadata"
TMP="$DB_DIR/tmp"
export KRAKEN2_DB="$DB_DIR"

mkdir -p "$DB_DIR/taxonomy" "$META" "$TMP" "$FA_DIR"

# Taxonomy (download once)
if [ ! -s "$DB_DIR/taxonomy/nodes.dmp" ]; then
  echo "[INFO] downloading taxonomy..."
  kraken2-build --download-taxonomy --skip-maps --db "$KRAKEN2_DB"
else
  echo "[INFO] taxonomy already present."
fi

# Check sanity
ls -lh "$DB_DIR/taxonomy"/{names.dmp,nodes.dmp,merged.dmp,delnodes.dmp}

# Build seqid→taxid map using accession lookups
MAP="$DB_DIR/taxonomy/seqid2taxid.map"
: > "$MAP"
echo -e "accession\tspecies\ttaxid" > "$META/species_taxids.tsv"

# Patch missing accession with distinct variable
MISSING_ACC="GCF_907164915.1"
if [ -n "$MISSING_ACC" ] && ! [ -s "$FA_DIR/${MISSING_ACC}.fna" ]; then
  echo "[INFO] fetching missing fasta for $MISSING_ACC"
  mkdir -p "$TMP/patch"
  datasets download genome accession "$MISSING_ACC" --include genome --filename "$TMP/patch/$MISSING_ACC.zip"
  unzip -q "$TMP/patch/$MISSING_ACC.zip" -d "$TMP/patch/$MISSING_ACC"
  cat "$TMP/patch/$MISSING_ACC"/ncbi_dataset/data/*/*.fna > "$FA_DIR/$MISSING_ACC.fna"
  TAXID_PATCH=$(datasets summary genome accession "$MISSING_ACC" | jq -r '.reports[0].organism.tax_id // empty')
  if [ -n "$TAXID_PATCH" ]; then
    awk -v T="$TAXID_PATCH" '/^>/{id=$1; sub(/^>/,"",id); print id "\t"T}' "$FA_DIR/$MISSING_ACC.fna" >> "$MAP"
    SPECIES_PATCH=$(datasets summary genome accession "$MISSING_ACC" | jq -r '.reports[0].organism.organism_name // empty')
    echo -e "$MISSING_ACC\t$SPECIES_PATCH\t$TAXID_PATCH" >> "$META/species_taxids.tsv"
  else
    echo "[WARN] no taxid for $MISSING_ACC"
  fi
  rm -rf "$TMP/patch"
fi

# Process the full accession list and append to map
while IFS= read -r ACC; do
  ACC="${ACC%%$'\r'}"
  [ -n "$ACC" ] || continue
  F="$FA_DIR/${ACC}.fna"
  if [ ! -s "$F" ]; then
    echo "[WARN] missing fasta for $ACC (skipping mapping)"; continue
  fi

  # Fetch taxid + species from accession metadata
  if datasets summary genome accession "$ACC" > "$TMP/${ACC}.json"; then
    TAXID=$(jq -r '.reports[0].organism.tax_id // empty' "$TMP/${ACC}.json")
    SPECIES=$(jq -r '.reports[0].organism.organism_name // empty' "$TMP/${ACC}.json")
  else
    TAXID=""; SPECIES=""
  fi
  if [ -z "$TAXID" ]; then
    echo "[WARN] no taxid for $ACC; skipping"
    continue
  fi
  echo -e "$ACC\t$SPECIES\t$TAXID" >> "$META/species_taxids.tsv"

# Append all sequence IDs -> taxid
  awk -v T="$TAXID" '/^>/{id=$1; sub(/^>/,"",id); print id "\t"T}' "$F" >> "$MAP"
done < "$ACC_LIST"

LINES=$(wc -l < "$MAP" || echo 0)
echo "[INFO] seqid2taxid.map rows: $LINES"
if [ "$LINES" -eq 0 ]; then
  echo "[ERROR] empty seqid2taxid.map; nothing to build."; exit 5
fi

# 3) Add library and build
echo "[INFO] adding FASTAs to Kraken2 library..."
shopt -s nullglob
for f in "$FA_DIR"/*.fna; do
  echo "[ADD] $(basename "$f")"
  kraken2-build --add-to-library "$f" --db "$KRAKEN2_DB"
done

echo "[INFO] building DB..."
kraken2-build --build --threads "$SLURM_CPUS_PER_TASK" --db "$KRAKEN2_DB"
kraken2-build --clean --db "$KRAKEN2_DB"

# 4) Sanity check
ls -lh "$DB_DIR"/{hash.k2d,opts.k2d,taxo.k2d}
kraken2-inspect --db "$KRAKEN2_DB" | head
echo "[DONE] DB at: $KRAKEN2_DB"
