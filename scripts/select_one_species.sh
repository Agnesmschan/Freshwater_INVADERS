#!/bin/bash
#SBATCH --partition=medium
#SBATCH --output=select_one_species.out
#SBATCH --error=select_one_species.err
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

rm -f ../data/refgen_one_each_species.tsv

awk 'BEGIN {
    OFS="\t";
    print "search", "species", "accession", "assembly_level";
}
NR > 1 {
    # assign assembly level priority (lower is better)
    if (NF < 4) next;
    level = tolower($4);
    if (level == "complete genome") rank = 1;
    else if (level == "chromosome") rank = 2;
    else if (level == "scaffold") rank = 3;
    else if (level == "contig") rank = 4;
    else rank = 99;

    key = $2;  # species name for uniqueness
    if (!(key in best) || rank < best_rank[key]) {
        best[key] = $0;
        best_rank[key] = rank;
    }
}
END {
    for (k in best) print best[k];
}' ../data/refgen_species.tsv > ../data/refgen_one_each_species.tsv


