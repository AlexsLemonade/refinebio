#!/bin/sh

# Thin wrapper around `docker buildx bake` for building a single image.
# Build configuration (tags, cache, platform, FROM resolution) lives in
# docker-bake.hcl; this script provides the equivalent CLI for callers that
# expect a per-image entry point.

set -e

script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# shellcheck disable=SC1091
. ./common.sh

# Run from the repo root so bake's relative context (".") resolves correctly.
cd ..

print_help() {
    cat <<EOF
Build a single docker image via 'docker buildx bake'.

Usage: $(basename "$0") -i IMAGE [-s SERVICE] [-r REPO] [-v VERSION] [-u] [-h]

Options:
  -i IMAGE     Image (bake target) to build, e.g. base, api, salmon. Required.
  -s SERVICE   Accepted for back-compat with existing callers; ignored.
               docker-bake.hcl already declares the dockerfile path per target.
  -r REPO      Override DOCKERHUB_REPO (default: ccdlstaging).
  -v VERSION   Override SYSTEM_VERSION (default: a hash of the current branch).
  -u           Push to DOCKERHUB_REPO instead of loading to the local daemon.
  -h           Print this help.
EOF
}

ACTION="--load"

while getopts "uhi:r:s:v:" opt; do
    case $opt in
    i) IMAGE_NAME="$OPTARG" ;;
    s) : ;; # accepted but no longer needed; see help
    r) export DOCKERHUB_REPO="$OPTARG" ;;
    v) export SYSTEM_VERSION="$OPTARG" ;;
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

if [ -z "$IMAGE_NAME" ]; then
    echo "Error: must specify image with -i" >&2
    exit 1
fi

# When pushing to a registry, also push cache layers (mode=max).
if [ "$ACTION" = "--push" ]; then
    export PUSH_CACHE=true
fi

docker buildx bake "$IMAGE_NAME" "$ACTION"

# When pushed, pull back so the image is also in the local daemon.
# Mirrors the original script's trailing `docker pull`.
if [ "$ACTION" = "--push" ]; then
    docker pull --platform linux/amd64 \
        "$DOCKERHUB_REPO/dr_$IMAGE_NAME:$SYSTEM_VERSION"
fi
