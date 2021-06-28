#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

while read -r row; do
    # Exporting an expansion rather than a variable, which is exactly what we want to do.
    # shellcheck disable=SC2163
    export "${row}"
done < ../infrastructure/prod_env

python3 kill_all_jobs.py
