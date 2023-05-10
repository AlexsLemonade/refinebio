#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

while read -r row; do
    # Exporting an expansion rather than a variable, which is exactly what we want to do.
    # shellcheck disable=SC2163
    export "${row}"
done <../infrastructure/prod_env

python3 kill_all_jobs.py
