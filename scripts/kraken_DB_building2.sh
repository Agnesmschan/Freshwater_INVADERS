#!/bin/bash
#SBATCH --job-name=kraken_DB_building2
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -euo pipefail

# (optional) ensure tools on PATH
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate Kraken2

# ---- paths ----
ACC_LIST=$(readlink -f ../data/genome_accessions_to_download.txt)
DB_DIR=$(readlink -f ../kraken2_databases/Freshwater_invaders)
FA_DIR="$DB_DIR/fastas"
META="$DB_DIR/metadata"
TMP="$DB_DIR/tmp"
export KRAKEN2_DB="$DB_DIR"

mkdir -p "$DB_DIR" "$FA_DIR" "$META" "$TMP" "$DB_DIR/taxonomy" "$DB_DIR/library"

# ---- taxonomy (minimal, no giant accession maps) ----
if [ ! -s "$DB_DIR/taxonomy/nodes.dmp" ]; then
  echo "[INFO] downloading taxonomy..."
  kraken2-build --download-taxonomy --skip-maps --db "$KRAKEN2_DB"
else
  echo "[INFO] taxonomy already present."
fi
ls -lh "$DB_DIR/taxonomy"/{names.dmp,nodes.dmp,merged.dmp,delnodes.dmp}

# ---- build seqid→taxid map (tab: seqid \t taxid), from FASTAs ----
MAP="$DB_DIR/taxonomy/seqid2taxid.map"
: > "$MAP"
echo -e "accession\tspecies\ttaxid" > "$META/species_taxids.tsv"

# Optional: patch a known-missing accession once.
MISSING_ACC="GCF_907164915.1"
if [ -n "$MISSING_ACC" ] && [ ! -s "$FA_DIR/${MISSING_ACC}.fna" ]; then
  echo "[INFO] fetching missing fasta for $MISSING_ACC"
  mkdir -p "$TMP/patch"
  datasets download genome accession "$MISSING_ACC" --include genome \
    --filename "$TMP/patch/$MISSING_ACC.zip"
  unzip -q "$TMP/patch/$MISSING_ACC.zip" -d "$TMP/patch/$MISSING_ACC"
  cat "$TMP/patch/$MISSING_ACC"/ncbi_dataset/data/*/*.fna > "$FA_DIR/$MISSING_ACC.fna"
  TAXID_PATCH=$(datasets summary genome accession "$MISSING_ACC" | jq -r '.reports[0].organism.tax_id // empty')
  SPECIES_PATCH=$(datasets summary genome accession "$MISSING_ACC" | jq -r '.reports[0].organism.organism_name // empty')
  if [ -n "$TAXID_PATCH" ]; then
    awk -v T="$TAXID_PATCH" '/^>/{id=$1; sub(/^>/,"",id); print id "\t"T}' \
      "$FA_DIR/$MISSING_ACC.fna" >> "$MAP"
    echo -e "$MISSING_ACC\t$SPECIES_PATCH\t$TAXID_PATCH" >> "$META/species_taxids.tsv"
  else
    echo "[WARN] no taxid for $MISSING_ACC"
  fi
  rm -rf "$TMP/patch"
fi

# Map each accession in the list (only if its FASTA exists)
while IFS= read -r ACC; do
  ACC="${ACC%%$'\r'}"; [ -n "$ACC" ] || continue
  F="$FA_DIR/${ACC}.fna"
  if [ ! -s "$F" ]; then
    echo "[WARN] missing fasta for $ACC (skipping mapping)"
    continue
  fi
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
  awk -v T="$TAXID" '/^>/{id=$1; sub(/^>/,"",id); print id "\t"T}' "$F" >> "$MAP"
done < "$ACC_LIST"

# Normalize line endings & expose the map where kraken2-build sees it
sed -i 's/\r$//' "$MAP"
ln -sf "$MAP" "$DB_DIR/seqid2taxid.map"
ln -sf "$MAP" "$DB_DIR/library/seqid2taxid.map"

# Remove any hints that force accession-map path
rm -f "$DB_DIR"/accmap_file.tmp "$DB_DIR"/seqid2taxid.map.tmp
rm -f "$DB_DIR/taxonomy"/accmap.dlflag "$DB_DIR/taxonomy"/nucl_*accession2taxid* 2>/dev/null || true

# Coverage check: every header should have a mapping
echo "[INFO] coverage check of map vs FASTA headers:"
ok=1
for f in "$FA_DIR"/*.fna; do
  n_hdr=$(grep -c '^>' "$f" || true)
  awk '/^>/{print substr($1,2)}' "$f" > "$TMP/ids.txt"
  n_map=$(grep -F -f "$TMP/ids.txt" "$MAP" | wc -l)
  printf "  %-20s headers=%-6s mapped=%-6s\n" "$(basename "$f")" "$n_hdr" "$n_map"
  if [ "$n_hdr" -ne "$n_map" ]; then ok=0; fi
done
if [ $ok -ne 1 ]; then
  echo "[ERROR] Some FASTA headers are not in seqid2taxid.map. Fix them then re-run."
  exit 6
fi

# ---- add FASTAs to library (one by one) ----
echo "[INFO] adding FASTAs to Kraken2 library..."
shopt -s nullglob
for f in "$FA_DIR"/*.fna; do
  echo "[ADD] $(basename "$f")"
  kraken2-build --add-to-library "$f" --db "$KRAKEN2_DB"
done

# ---- build & clean ----
threads="${SLURM_CPUS_PER_TASK:-4}"
echo "[INFO] building DB with $threads threads..."
kraken2-build --build --threads "$threads" --db "$KRAKEN2_DB"
kraken2-build --clean --db "$KRAKEN2_DB"

# ---- sanity ----
ls -lh "$DB_DIR"/{hash.k2d,opts.k2d,taxo.k2d} || true
kraken2-inspect --db "$KRAKEN2_DB" | head || true
echo "[DONE] DB at: $KRAKEN2_DB"
