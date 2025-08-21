#!/bin/bash
#SBATCH --partition=medium
#SBATCH --output=get_accessions.out
#SBATCH --error=get_accessions.err
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

mkdir -p ../data/genome_signatures
cat > ../data/genome_signatures/place_holder.sig

tail -n +2 ../data/refgen_one_each_species.tsv | awk -F'\t' '$3 != "None" { print $3 }' > ../data/genome_accessions.txt

comm -23 <(sort ../data/genome_accessions.txt) <(ls ../data/genome_signatures/*.sig | sed 's#.*/##; s/\.sig//g' | sort) > ../data/genome_accessions_to_download.txt

rm -f  ../data/genome_signatures/place_holder.sig
