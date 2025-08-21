#!/bin/bash

#SBATCH --output=rename_sig.out
#SBATCH --error=rename_sig.err
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Input files
MAPPING_FILE="../data/refgen_one_each_species.tsv"
SIG_DIR="../data/genome_signatures"
OUT_DIR="../data/Invaders_genome_signatures"

# Create output directory
mkdir -p "$OUT_DIR"

# Skip header line and read TSV
tail -n +2 "$MAPPING_FILE" | while IFS=$'\t' read -r search species accession assembly_level; do
    sig_path="${SIG_DIR}/${accession}_k31.sig"

    # Clean species name: replace spaces with underscores
    species_clean=$(echo "$species" | tr ' ' '_')

    out_sig="${OUT_DIR}/${accession}_${species_clean}_k31.sig"

    # Check if the .sig file exists
    if [[ -f "$sig_path" ]]; then
        echo "Renaming: $accession --> $species"
        sourmash sig rename "$sig_path" "$species" -o "$out_sig"
    else
        echo "WARNING: Signature not found for accession $accession"
    fi
done
