#!/bin/sh -e

# This script downloads the json file published by MetaSRA and translates it
# into the format that our metadata processing code understands.

METASRA_URL="https://metasra.biostat.wisc.edu/static/metasra_versions/v1.8/metasra.v1-8.json"
JSON_FILE="$(basename "$METASRA_URL")"
ACCESSIONS_URL="https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
ACCESSIONS_FILE="$(basename "$ACCESSIONS_URL")"
PROCESSED_ACCESSIONS_FILE="$(basename -s .tab "$ACCESSIONS_FILE").processed.tab"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

if ! [ -e "$JSON_FILE" ]; then
    wget --no-check-certificate "$METASRA_URL"
fi

if ! [ -e "$PROCESSED_ACCESSIONS_FILE" ]; then
    [ -e "$ACCESSIONS_FILE" ] || wget "$ACCESSIONS_URL"

    # Do some pre-processing on the accessions file so that it can fit in memory
    # in the python program.
    awk -f filter_accessions.awk "$ACCESSIONS_FILE" > "$PROCESSED_ACCESSIONS_FILE"

    # Clean up this big file once we've extracted the information we need
    rm "$PROCESSED_ACCESSIONS_FILE"
fi

python3 translate.py "$JSON_FILE" "$PROCESSED_ACCESSIONS_FILE"
