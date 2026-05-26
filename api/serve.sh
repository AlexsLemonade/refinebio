#!/bin/sh

# Serve the API on localhost:8000.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Run from the repo root so bin/rbio and compose.yml resolve.
cd ..

./bin/rbio build api_local

exec docker compose run --rm --service-ports api \
    python3 manage.py runserver 0.0.0.0:8000 "$@"
