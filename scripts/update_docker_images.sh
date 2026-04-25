#!/bin/sh

# Thin wrapper around `docker buildx bake` for building all images at once.
# Build configuration lives in docker-bake.hcl. Also (re)builds the common
# sdist that worker/foreman dockerfiles COPY in.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# shellcheck disable=SC1091
. ./common.sh

# Run from the repo root so bake's relative context and sdist paths resolve.
cd ..

print_help() {
    cat <<EOF
Build all docker images via 'docker buildx bake'.

Usage: $(basename "$0") [-r REPO] [-v VERSION] [-a] [-u] [-h]

Options:
  -r REPO      Override DOCKERHUB_REPO (default: ccdlstaging).
  -v VERSION   Override SYSTEM_VERSION (default: a hash of the current branch).
  -a           Also build affymetrix (uses bake's "all" group; default omits it).
  -u           Push to DOCKERHUB_REPO instead of loading to the local daemon.
  -h           Print this help.
EOF
}

ACTION="--load"
GROUP="default"

while getopts ":r:v:auh" opt; do
    case $opt in
    r) export DOCKERHUB_REPO="$OPTARG" ;;
    v) export SYSTEM_VERSION="$OPTARG" ;;
    a) GROUP="all" ;;
    u) ACTION="--push" ;;
    h)
        print_help
        exit 0
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        print_help >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
    esac
done

if [ "$ACTION" = "--push" ]; then
    export PUSH_CACHE=true
fi

# Set the version for the common project, then build the sdist that
# worker/foreman dockerfiles COPY in. Both are preconditions for those builds.
echo "$SYSTEM_VERSION" >common/version
rm -f common/dist/*
(cd common && python3 setup.py sdist 1>/dev/null)

docker buildx bake "$GROUP" "$ACTION"
