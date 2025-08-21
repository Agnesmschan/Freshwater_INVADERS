#!/bin/bash

#SBATCH --partition=medium
#SBATCH --job-name=check_target_species
#SBATCH --output=check_target_species.out
#SBATCH --error=check_target_species.err
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

# load ncbi tools
source ~/miniforge3/etc/profile.d/conda.sh
conda activate smash

# Input and output files
INPUT_SPECIES="../data/Target_species.txt"
OUT_DIR="../data"
LOG_DIR="logs"
TEMP_JSON="temp.json"

mkdir -p "$OUT_DIR" "$LOG_DIR"


# Check for existing reference genomes

echo -e "search\tspecies\taccession\tassembly_level" > "$OUT_DIR/refgen_species.tsv"

while IFS= read -r species; do
      echo "[INFO] Searching for genome of: $species"
      datasets summary genome taxon "$species" > "$TEMP_JSON" || true

      count=$(jq '.reports | length' "$TEMP_JSON")
      if [[ "$count" -gt 0 ]]; then
          jq -r --arg search "$species" '
              .reports[]? |
              [
                   $search,
                   .organism.organism_name,
                   .accession,
                   .assembly_info.assembly_level
              ] | @tsv
          ' "$TEMP_JSON" >> "$OUT_DIR/refgen_species.tsv"
      else
          echo -e "$species\t$species\tNone\tNone" >> "$OUT_DIR/refgen_species.tsv"
      fi
      rm -f "$TEMP_JSON"
done < "$INPUT_SPECIES"


# filtering results


awk -F'\t' 'NR > 1 && $3 == "None" { print $1}' "$OUT_DIR/refgen_species.tsv" > "$OUT_DIR/species_without_refgens.txt"
awk -F'\t' 'NR > 1 && $3 != "None" { print $1}' "$OUT_DIR/refgen_species.tsv" | sort | uniq > "$OUT_DIR/species_with_refgens.txt"
awk -F'\t' 'NR == 1 || $1 != $2' "$OUT_DIR/refgen_species.tsv" > "$OUT_DIR/mismatched_species.txt"

# Get family name for species with genomes

> "$OUT_DIR/species_family.txt"

while IFS= read -r species; do
    echo "[INFO] Getting family for: $species"
    datasets summary taxonomy taxon "$species" > "$TEMP_JSON" || true

    family=$(jq -r '.reports[0].taxonomy.classification.family.name // "None"' "$TEMP_JSON")
    echo "$family" >> "$OUT_DIR/species_family.txt"

    rm -f "$TEMP_JSON"
done < "$OUT_DIR/species_with_refgens.txt"

sort "$OUT_DIR/species_family.txt" | uniq > "$OUT_DIR/unique_species_family.txt"


# Search for genomes in those families

echo -e "search\tspecies\taccession\tassembly_level" > "$OUT_DIR/refgen_family.tsv"

while IFS= read -r family; do
    echo "[INFO] Searching for genomes in family: $family"
    datasets summary genome taxon "$family" > "$TEMP_JSON" || true

    count=$(jq '.reports | length' "$TEMP_JSON")
    if [[ "$count" -gt 0 ]]; then
        jq -r --arg search "$family" '
            .reports[]? |
            [
                $search,
                .organism.organism_name,
                .accession,
                .assembly_info.assembly_level
            ] | @tsv
         ' "$TEMP_JSON" >> "$OUT_DIR/refgen_family.tsv"
       else
           echo -e "$family\t$family\tNone\tNone" >> "$OUT_DIR/refgen_family.tsv"
       fi

       rm -f "$TEMP_JSON"
done < "$OUT_DIR/unique_species_family.txt"
