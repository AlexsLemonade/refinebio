#!/bin/sh -e

# This script downloads the json file published by MetaSRA and translates it
# into the format that our metadata processing code understands.

METASRA_URL="http://metasra.biostat.wisc.edu/static/metasra_versions/v1.6/metasra.v1-6.json"
JSON_FILE="$(basename "$METASRA_URL")"
SRADB_FILE="SRAmetadb.sqlite"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

if ! [ -e "$JSON_FILE" ]; then
    wget "$METASRA_URL"
fi

# We are downlading the 24GB SRAmetadb.sqlite to use it to translate SRS's to
# SRR's. This seems overkill, but the alternative is hammering ENA for every
# single accession in the MetaSRA file which could take *at least* a couple
# hours.
if ! [ -e "$SRADB_FILE" ]; then
    wget "https://s3.amazonaws.com/starbuck1/sradb/SRAmetadb.sqlite.gz"
    gunzip "$SRADB_FILE.gz"
fi

python3 translate.py "$JSON_FILE"
