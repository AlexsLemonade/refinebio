#!/bin/sh

# Thin shim around `docker compose up -d postgres`.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Run from the repo root so compose.yml resolves.
cd ..

exec docker compose up -d postgres
