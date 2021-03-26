#!/bin/sh

# Runs all common and app specific tests.

# Don't keep going if some of the tests fail
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of the refinebio project
cd ..

# First ensure postgres is running
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
# Then ensure elasticsearch is running
elif ! [ "$(docker ps --filter name=dres -q)" ]; then
    echo "You must start Elasticsearch first with:" >&2
    echo "./scripts/run_es.sh" >&2
    exit 1
fi

mkdir -p test_volume && chmod -R a+rw test_volume
./scripts/update_models.sh
./api/run_tests.sh
./common/run_tests.sh
./foreman/run_tests.sh
./workers/run_tests.sh
