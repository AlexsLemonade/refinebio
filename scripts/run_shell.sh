#!/bin/sh

# Open an interactive Python shell in the foreman container.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Set up the foreman data volume directory if it does not already exist.
volume_directory="$script_directory/../foreman/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

# Run from the repo root so prepare_image.sh and compose.yml resolve.
cd ..

./scripts/prepare_image.sh -i foreman -s foreman

exec docker compose run --rm foreman python3 manage.py shell
