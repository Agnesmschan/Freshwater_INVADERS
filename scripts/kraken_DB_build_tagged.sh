#!/bin/bash
#SBATCH --job-name=kraken_DB_build_tagged_v3
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -euo pipefail

# --- paths ---
ACC_LIST=$(readlink -f ../data/genome_accessions_to_download.txt)  # one accession per line
DB_DIR=$(readlink -f ../kraken2_databases/Freshwater_invaders)
FA_DIR="$DB_DIR/fastas"                 # untagged FASTAs
TAG_DIR="$DB_DIR/fastas_tagged"         # tagged FASTAs
META="$DB_DIR/metadata"
TMP="$DB_DIR/tmp"
export KRAKEN2_DB="$DB_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$DB_DIR" "$FA_DIR" "$TAG_DIR" "$META" "$TMP" "$DB_DIR/library" "$DB_DIR/taxonomy"

# --- taxonomy (no accession maps) ---
if [ ! -s "$DB_DIR/taxonomy/nodes.dmp" ]; then
  echo "[INFO] downloading taxonomy (no accession maps)..."
  kraken2-build --download-taxonomy --skip-maps --db "$KRAKEN2_DB"
else
  echo "[INFO] taxonomy already present."
fi

# --- tag FASTA headers with |kraken:taxid|T ---
echo -e "accession\tspecies\ttaxid" > "$META/species_taxids.tsv"
shopt -s nullglob
for F in "$FA_DIR"/*.fna; do
  ACC=$(basename "$F" .fna)
  if ! datasets summary genome accession "$ACC" > "$TMP/${ACC}.json"; then
    echo "[WARN] metadata fetch failed for $ACC; skipping"; continue
  fi
  TAXID=$(jq -r '.reports[0].organism.tax_id // empty' "$TMP/${ACC}.json")
  SPECIES=$(jq -r '.reports[0].organism.organism_name // empty' "$TMP/${ACC}.json")
  if [ -z "$TAXID" ]; then echo "[WARN] no taxid for $ACC; skipping"; continue; fi
  echo -e "$ACC\t$SPECIES\t$TAXID" >> "$META/species_taxids.tsv"

  TGT="$TAG_DIR/${ACC}.tag.fna"
  awk -v T="$TAXID" '
    /^>/ {
      sub(/\r$/, "");
      if (index($0, "kraken:taxid|") == 0) print $0"|kraken:taxid|"T;
      else print $0;
      next
    } { print }
  ' "$F" > "$TGT"
done

# --- build seqid→taxid.map FROM TAGGED FASTAS; place where Kraken expects ---
MAP="$DB_DIR/taxonomy/seqid2taxid.map"
: > "$MAP"
for F in "$TAG_DIR"/*.tag.fna; do
  awk '
    /^>/ {
      tid="";
      if (match($0, /\|kraken:taxid\|([0-9]+)/, m)) tid=m[1];
      id=$1; sub(/^>/,"",id);
      if (tid != "") print id "\t" tid;
    }
  ' "$F" >> "$MAP"
done
sed -i 's/\r$//' "$MAP"
# COPY (not symlink) to exact locations the build script tests/uses:
cp -f "$MAP" "$DB_DIR/seqid2taxid.map"
cp -f "$MAP" "$DB_DIR/library/seqid2taxid.map"

# Debug: show where map is and size
echo "[DEBUG] maps:"
ls -lh "$DB_DIR"/seqid2taxid.map "$DB_DIR/library/seqid2taxid.map" "$DB_DIR/taxonomy/seqid2taxid.map"

# --- remove anything that causes accession-map path or build-skip ---
rm -f "$DB_DIR"/accmap_file.tmp "$DB_DIR"/seqid2taxid.map.tmp
rm -f "$DB_DIR/taxonomy"/accmap.dlflag "$DB_DIR/taxonomy"/nucl_*accession2taxid* 2>/dev/null || true
rm -f "$DB_DIR"/hash.k2d "$DB_DIR"/opts.k2d "$DB_DIR"/taxo.k2d

# --- add TAGGED FASTAs with NO MASKING (don’t alter headers) ---
echo "[INFO] adding tagged FASTAs to library..."
for F in "$TAG_DIR"/*.tag.fna; do
  echo "[ADD] $(basename "$F")"
  kraken2-build --add-to-library "$F" --no-masking --db "$KRAKEN2_DB"
done

# Confirm library has files and prelim_map.txt exists (good sign)
echo "[DEBUG] library contents:"
find "$DB_DIR/library" -maxdepth 2 -type f | head
echo "[DEBUG] prelim_map.txt (if any):"
find "$DB_DIR/library" -maxdepth 2 -name 'prelim_map.txt' -exec head -2 {} \; -print || true

# --- build & clean ---
echo "[INFO] building DB with $THREADS threads..."
kraken2-build --build --threads "$THREADS" --db "$KRAKEN2_DB"
kraken2-build --clean --db "$KRAKEN2_DB"

echo "[INFO] DB files:"
ls -lh "$DB_DIR"/{hash.k2d,opts.k2d,taxo.k2d} || true
echo "[INFO] kraken2-inspect (head):"
kraken2-inspect --db "$KRAKEN2_DB" | head || true

