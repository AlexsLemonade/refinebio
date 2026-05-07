#!/bin/sh

# Serve the API on localhost:8000.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Run from the repo root so prepare_image.sh and compose.yml resolve.
cd ..

./scripts/prepare_image.sh -i api_local -s api

exec docker compose run --rm --service-ports api \
    python3 manage.py runserver 0.0.0.0:8000 "$@"
