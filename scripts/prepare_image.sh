#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# shellcheck disable=SC1091
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
    DOCKER_ACTION="--load"
fi

DOCKER_FILE_PATH="$SERVICE/dockerfiles/Dockerfile.$IMAGE_NAME"

attempt=0
max_attempts=3
finished=1
while [ $finished != 0 ] && [ $attempt -lt $max_attempts ]; do
    if [ $attempt -gt 0 ]; then
        echo "Failed to build $IMAGE_NAME:$SYSTEM_VERSION image, trying again."
    fi

    echo "Building the $IMAGE_NAME:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
    update_docker_image "$DOCKERHUB_REPO" "$IMAGE_NAME" "$SYSTEM_VERSION" "$DOCKER_FILE_PATH" "$DOCKER_ACTION"

    finished=$?
    attempt=$((attempt + 1))
done

if [ $finished -ne 0 ] && [ $attempt -ge $max_attempts ]; then
    echo "Could not build $DOCKERHUB_IMAGE after $attempt attempts."
    exit 1
fi

if test "$GITHUB_ACTION"; then
    docker pull \
        --platform linux/amd64 \
        "$DOCKERHUB_IMAGE"
fi
