#!/bin/bash

#SBATCH --output=compare_k_mer.out
#SBATCH --error=compare_k_mer.err
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1

set -euo pipefail

# Compare k=31
echo "Comparing k=31 signatures..."
sourmash compare ../data/genome_signatures/*_k31.sig -o k31_cmp

echo "Plotting k=31..."
sourmash plot k31_cmp --labels --output-dir ../data/k31_plots

# Compare k=51
echo "Comparing k=51 signatures..."
sourmash compare ../data/genome_signatures/*_k51.sig -o k51_cmp

echo "Plotting k=51..."
sourmash plot k51_cmp --labels --output-dir ../data/k51_plots

