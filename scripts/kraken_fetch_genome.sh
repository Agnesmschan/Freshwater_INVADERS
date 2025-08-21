#!/bin/bash

#SBATCH --job-name=kraken_fetch_genome
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=short

source ~/miniforge3/etc/profile.d/conda.sh
conda activate Kraken2

set -euo pipefail

TSV=$(readlink -f ../data/refgen_one_each_species.tsv)
DB_DIR=$(readlink -f ../kraken2_databases/Freshwater_invaders)
OUT_FASTA="$DB_DIR/fastas"
OUT_META="$DB_DIR/metadata"
TMP="$DB_DIR/tmp"

export KRAKEN2_DB="$DB_DIR"

mkdir -p "$DB_DIR" "$OUT_FASTA" "$OUT_META" "$TMP" "$DB_DIR/taxonomy"

# Normalise TSV newlines to avoid CRLF surprises (writes a temp copy)
TSV_NORM="$TMP/targets.tsv"
sed 's/\r$//' "$TSV" > "$TSV_NORM"

# Set map and download taxonomy once
if [ ! -s "$DB_DIR/taxonomy/nodes.dmp" ]; then
  kraken2-build --download-taxonomy --db "$KRAKEN2_DB"
fi

MAP="$DB_DIR/taxonomy/seqid2taxid.map"
: > "$MAP"
echo -e "accession\tspecies\ttaxid" > "$OUT_META/species_taxids.tsv"


# Read tsv file, skip header and handle None fields
tail -n +2 "$TSV_NORM" | awk -F'\t' 'NF>=3 {print $2 "\t" $3}' | \
while IFS=$'\t' read -r species acc; do
  if [[ -z "${acc:-}" || "$acc" == "None" ]]; then
    echo "[WARN] skipping '$species' (no accession)"
    continue
  fi

  echo "[INFO] $acc  $species"

# Get taxid from accession
taxid=""
  if datasets summary genome accession "$acc" > "$TMP/${acc}.json"; then
    taxid=$(jq -r '.reports[0].organism.tax_id // empty' "$TMP/${acc}.json")
  fi
  if [[ -z "$taxid" || "$taxid" == "null" ]]; then
    if datasets summary taxonomy taxon "$species" > "$TMP/${acc}.tax.json"; then
      taxid=$(jq -r '.reports[0].tax_id // empty' "$TMP/${acc}.tax.json")
    fi
  fi
  if [[ -z "$taxid" ]]; then
    echo "[WARN] no taxid for $acc ($species); skipping"
    continue
  fi
  echo -e "$acc\t$species\t$taxid" >> "$OUT_META/species_taxids.tsv"

# Download genome FASTA(s) for that accession
  if ! datasets download genome accession "$acc" --include genome --filename "$TMP/${acc}.zip"; then
    echo "[WARN] download failed for $acc"
    continue
  fi
  unzip -q "$TMP/${acc}.zip" -d "$TMP/${acc}/"

# Concatenate all .fna into one file per accession
  out_fa="$OUT_FASTA/${acc}.fna"
  find "$TMP/${acc}/ncbi_dataset/data" -type f -name '*.fna' -exec cat {} + > "$out_fa" || true
  if [ ! -s "$out_fa" ]; then
    echo "[WARN] no .fna for $acc; skipping"
    rm -rf "$TMP/${acc}" "$TMP/${acc}.zip"
    continue
  fi

# Append seqid→taxid mappings for all contigs
  awk -v T="$taxid" '/^>/{id=$1; sub(/^>/,"",id); print id "\t"T}' "$out_fa" >> "$MAP"

# Add to Kraken2 library
kraken2-build --add-to-library "$out_fa" --db "$KRAKEN2_DB"

# Cleanup tmp for this accession
  rm -rf "$TMP/${acc}" "$TMP/${acc}.zip"
done

# Build the DB
echo "[INFO] building Kraken2 DB..."
kraken2-build --build --threads "$SLURM_CPUS_PER_TASK" --db "$KRAKEN2_DB"
kraken2-build --clean --db "$KRAKEN2_DB"

# Finally... Sanity check
echo "[INFO] verifying DB..."
ls -lh "$DB_DIR"/{hash.k2d,opts.k2d,taxo.k2d}
kraken2-inspect --db "$KRAKEN2_DB" | head
echo "[DONE] DB at: $KRAKEN2_DB"
echo "[INFO] FASTAs: $OUT_FASTA"
echo "[INFO] Taxid map: $MAP"
