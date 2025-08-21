#!/bin/bash

#SBATCH --partition=short
#SBATCH --output=get_sra_accessions.out
#SBATCH --error=get_sra_accessions.err
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

set -euo pipefail

# Set 3 user input, modify for custom searches:

SEARCH_TEXT="((((United Kingdom[All Fields] OR UK[Title]) AND eDNA[All Fields]) AND (Environmental[All Fields] AND DNA[All Fields])) AND Freshwater[All Fields])"
SEARCH_NAME="UK_freshwater_eDNA"
UPDATE=false #true or false

# set directory
DATA_DIR="../Freshwater/data"
OUTFILE="${DATA_DIR}/${SEARCH_NAME}.txt"
mkdir -p "$DATA_DIR"


# Search SRA and find accessions
if [ "$UPDATE" = true ]; then
    echo "Updating previous results for '${SEARCH_NAME}'..."

    if [ ! -f "$OUTFILE" ]; then
       echo "[ERROR] NO existing file to update. Set UPDATE=false or check filename."
       exit 1
    fi

COUNT_BEFORE=$(wc -l < "$OUTFILE")
DATE_VAR=$(stat -c "%y" "$OUTFILE" | cut -d ' ' -f1)
BACKUP_FILE="${DATA_DIR}/${SEARCH_NAME}_before_${DATE_VAR}.txt"

cp "$OUTFILE" "$BACKUP_FILE"
> "$OUTFILE"

esearch -db sra -query "$SEARCH_TEXT" | \
  efilter -mindate "$DATE_VAR" | \
  efetch -format runinfo | \
  cut -d "," -f 1 | sed '1d' >> "$OUTFILE"

COUNT_AFTER=$(wc -l < "$OUTFILE")
   echo "Updated $SEARCH_NAME: $COUNT_BEFORE to $COUNT_AFTER accessions."

else

   echo "Running new search for '${SEARCH_NAME}'..."
   > "$OUTFILE"

esearch -db sra -query "$SEARCH_TEXT" | \
  efetch -format runinfo | \
  cut -d "," -f 1 | sed '1d' >> "$OUTFILE"

COUNT=$(wc -l < "$OUTFILE")
   echo "Found ${COUNT} accessions for query: ${SEARCH_TEXT}"

fi
