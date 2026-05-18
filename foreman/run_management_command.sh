#!/bin/sh

# Run a Django management command in the foreman container.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Set up the data volume directory if it does not already exist.
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

# Run from the repo root so compose.yml resolves.
cd ..

./bin/rbio build --load foreman

exec docker compose run --rm foreman python3 manage.py "$@"
