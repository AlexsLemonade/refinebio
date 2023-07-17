#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

. ./common.sh

# We need access to all of the projects
cd ..

print_description() {
    echo "Prepares an image specified by -i."
    echo "Must be called from the repo root."
}

print_options() {
    echo "Options:"
    echo "    -h           Prints the help message."
    echo "    -i IMAGE     The image to be prepared. This must be specified."
    echo "    -s SERVICE   The service to seach for a dockerfile."
    echo "                 The default option is 'workers'."
    echo "    -u           Push the built image to the Dockerhub."
    echo "    -r REPO      The docker registry to use for pull/push actions."
    echo "                 The default option is 'ccdlstaging'."
    echo
    echo "Examples:"
    echo "    Build the image ccdl/dr_downloaders:"
    echo "    ./scripts/prepare_image.sh -i downloaders -r ccdlstaging"
}

while getopts "uhi:r:s:" opt; do
    case $opt in
    i)
        IMAGE_NAME="$OPTARG"
        ;;
    r)
        DOCKERHUB_REPO="$OPTARG"
        ;;

    s)
        SERVICE="$OPTARG"
        ;;
    u)
        DOCKER_ACTION="--push"
        ;;
    h)
        print_description
        echo
        print_options
        exit 0
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        print_options >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        print_options >&2
        exit 1
        ;;
    esac
done

if [ -z "$IMAGE_NAME" ]; then
    echo "Error: you must specify an image with -i" >&2
    exit 1
fi

if [ -z "$SERVICE" ]; then
    SERVICE="workers"
fi

if [ -z "$DOCKERHUB_REPO" ]; then
    DOCKERHUB_REPO="ccdlstaging"
fi

# Defaults to commit hash value for if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="$(get_branch_hash)"
fi

if [ -z "$DOCKER_ACTION" ]; then
    DOCKER_ACTION="--push"
fi

DOCKERHUB_IMAGE="$DOCKERHUB_REPO/dr_$IMAGE_NAME"

CACHE_FROM_LATEST="type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest"
CACHE_FROM_VERSION="type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION"
CACHE_TO_LATEST="type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest,mode=max"
CACHE_TO_VERSION="type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION,mode=max"

if test "$GITHUB_ACTION"; then
    CACHE_TO_LATEST="type=gha"
    CACHE_TO_VERSION="type=gha"
    DOCKER_ACTION="--push"
fi

DOCKER_FILE_PATH="$SERVICE/dockerfiles/Dockerfile.$IMAGE_NAME"

echo
echo "Building the $IMAGE_NAME:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
echo

attempt=0
attempts=3
finished=1
while [ $finished != 0 ] && [ $attempt -lt $attempts ]; do
    if [ $attempt -gt 0 ]; then
        echo "Failed to build $IMAGE_NAME:$SYSTEM_VERSION image, trying again."
    fi

    set_up_docker_builder

    docker buildx build \
        --build-arg DOCKERHUB_REPO="$DOCKERHUB_REPO" \
        --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
        --cache-from "$CACHE_FROM_LATEST" \
        --cache-from "$CACHE_FROM_VERSION" \
        --cache-to "$CACHE_TO_LATEST" \
        --cache-to "$CACHE_TO_VERSION" \
        --file "$DOCKER_FILE_PATH" \
        --platform linux/amd64 \
        --tag "$DOCKERHUB_IMAGE:latest" \
        --tag "$DOCKERHUB_IMAGE:$SYSTEM_VERSION" \
        "$DOCKER_ACTION" \
        .

    finished=$?
    attempt=$((attempt + 1))
done

if [ $finished -ne 0 ] && [ $attempt -ge $attempts ]; then
    echo "Could not build $DOCKERHUB_IMAGE after $attempt attempts."
    exit 1
fi

docker pull \
    --platform linux/amd64 \
    "$DOCKERHUB_IMAGE:$SYSTEM_VERSION"
