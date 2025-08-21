#!/bin/bash
#SBATCH --job-name=unpack_core_nt
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=unpack_core_nt.out
#SBATCH --error=unpack_core_nt.err

set -euo pipefail
DEST="/mnt/shared/projects/nhm/clark-student/AirSeq/Reference_Databases/core_nt_20250609"
cd "$DEST"

mkdir -p core_nt_database
# use pigz if available for speed
if command -v pigz >/dev/null 2>&1; then
  tar -I pigz -xf k2_core_nt_20250609.tar.gz -C core_nt_20250609
else
  tar -xzf k2_core_nt_20250609.tar.gz -C core_nt_20250609 --checkpoint=. --totals
fi

# sanity checks
ls -lh core_nt_20250609/{hash.k2d,opts.k2d,taxo.k2d}
ls -lh core_nt_20250609/taxonomy/{names.dmp,nodes.dmp,merged.dmp,delnodes.dmp} || true
kraken2-inspect --db core_nt_20250609 | head
EOF

sbatch unpack_core_nt.sbatch
