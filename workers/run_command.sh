#!/bin/sh

# Run a command inside a worker container.

set -e

while getopts "i:" opt; do
    case $opt in
    i)
        IMAGE="$OPTARG"
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    esac
done

if [ -z "$IMAGE" ]; then
    IMAGE="smasher"
else
    shift
    shift
fi

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Ensure that postgres is running.
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

# Run from the repo root so prepare_image.sh and compose.yml resolve.
cd ..

./scripts/prepare_image.sh -i "$IMAGE" -s workers

exec docker compose run --rm "$IMAGE" bash -c "$@"
