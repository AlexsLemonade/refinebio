#!/bin/sh

# Runs all common and app specific tests.

# Don't keep going if some of the tests fail
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Get access to all of the refinebio project
cd ..

# First ensure postgres is running
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./bin/rbio dev:up -s postgres" >&2
    exit 1
# Then ensure elasticsearch is running
elif ! [ "$(docker ps --filter name=dres -q)" ]; then
    echo "You must start Elasticsearch first with:" >&2
    echo "./bin/rbio dev:up -s elasticsearch" >&2
    exit 1
fi

mkdir -p test_volume && chmod -R a+rw test_volume
./bin/rbio common:update
./bin/rbio test:api
./bin/rbio test:common
./bin/rbio test:foreman
./workers/run_tests.sh
