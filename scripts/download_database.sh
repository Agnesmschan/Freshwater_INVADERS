#!/bin/bash
#SBATCH --job-name=downlaod_core_nt
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=download_core_nt.out
#SBATCH --error=download_core_nt.err

set -euo pipefail

URL="https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20250609.tar.gz"
MD5="https://genome-idx.s3.amazonaws.com/kraken/core_nt_20250609/core_nt.md5"
DEST="/mnt/shared/projects/nhm/clark-student/AirSeq/Reference_Databases/core_nt_20250609"

mkdir -p "$DEST"
cd "$DEST"

echo "[INFO] free space:"
df -h .

# robust/resumable download (wget)
  wget -c --show-progress "$URL"
  wget -c --show-progress "$MD5"

echo "[INFO] verifying checksum..."
md5sum -c k2_core_nt_20250609.md5

echo "[OK] download complete and verified."
