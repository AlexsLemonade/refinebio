#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Import the functions in common.sh
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
    echo "    -d           Pull the latest version of the image from Dockerhub."
    echo "    -u           Push the built image to the Dockerhub."
    echo "    -r REPO      The docker registry to use for pull/push actions."
    echo "                 The default option is 'ccdlstaging'."
    echo
    echo "Examples:"
    echo "    Build the image ccdl/dr_downloaders:"
    echo "    ./scripts/prepare_image.sh -i downloaders -r ccdlstaging"
}

while getopts "udhi:b:r:s:" opt; do
    case $opt in
    b)
        DOCKER_BUILDER="$OPTARG"
        ;;
    d)
        PULL="True"
        ;;
    i)
        IMAGE="$OPTARG"
        ;;
    r)
        DOCKERHUB_REPO="$OPTARG"
        ;;

    s)
        SERVICE="$OPTARG"
        ;;
    u)
        DOCKER_ACTION="push"
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

if [ -z "$IMAGE" ]; then
    echo "Error: you must specify an image with -i" >&2
    exit 1
fi

if [ -z "$SERVICE" ]; then
    SERVICE="workers"
fi

if [ -z "$DOCKERHUB_REPO" ]; then
    DOCKERHUB_REPO="ccdlstaging"
fi

# Default to "local" for system version if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="local$(date +%s)"
fi

if [ -z "$DOCKER_ACTION" ]; then
    DOCKER_ACTION="load"
fi

# We want to check if a test image has been built for this branch. If
# it has we should use that rather than building it slowly.
DOCKER_IMAGE="$DOCKERHUB_REPO/dr_$IMAGE"
# shellcheck disable=SC2086
if [ "$(docker_img_exists $IMAGE_NAME $branch_name)" ]; then
    docker pull "$IMAGE_NAME:$branch_name"
elif [ -n "$PULL" ]; then
    docker pull \
        --platform linux/amd64 \
        "$DOCKER_IMAGE"
else
    echo ""
    echo "Rebuilding the $DOCKER_IMAGE image."

    attempts=0
    finished=1
    while [ $finished != 0 ] && [ $attempts -lt 3 ]; do
        if [ $attempts -gt 0 ]; then
            echo "Failed to build $DOCKER_IMAGE, trying again."
        fi

        CURRENT_IMAGE="$DOCKER_IMAGE:$SYSTEM_VERSION"
        LATEST_IMAGE="$DOCKER_IMAGE:latest"
        if test "$GITHUB_ACTIONS"; then
            # Docker needs repositories to be lowercase.
            CACHE_REPO="$(echo "ghrc.io/$GITHUB_REPOSITORY" |
                tr '[:upper:]' '[:lower:]')"
            LATEST_IMAGE="$CACHE_REPO/dr_$IMAGE"
        fi

        if test "$DOCKER_BUILDER"; then
            docker buildx use "$DOCKER_BUILDER"
        fi

        DOCKER_BUILDKIT=1 docker buildx build \
            --"$DOCKER_ACTION" \
            --build-arg BUILDKIT_INLINE_CACHE=1 \
            --build-arg DOCKERHUB_REPO="$DOCKERHUB_REPO" \
            --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
            --cache-from "$LATEST_IMAGE" \
            --cache-from="$CURRENT_IMAGE" \
            --file "$SERVICE/dockerfiles/Dockerfile.$IMAGE" \
            --platform linux/amd64 \
            --tag "$DOCKER_IMAGE" \
            .
        finished=$?
        attempts=$((attempts + 1))
    done

    if [ $finished != 0 ] && [ $attempts -ge 3 ]; then
        echo "Could not build $DOCKER_IMAGE after three attempts."
        exit 1
    fi
fi
