#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
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
    echo "    -h           Prints the help message"
    echo "    -i IMAGE     The image to be prepared. This must be specified."
    echo "    -s SERVICE   The service to seach for a dockerfile."
    echo "                 The default option is 'workers'"
    echo "    -p           Pull the latest version of the image from Dockerhub"
    echo "    -d REPO      The docker repo to pull images from."
    echo "                 The default option is 'ccdl'"
    echo "    -b           The remote docker builder name (optional)."
    echo
    echo "Examples:"
    echo "    Build the image ccdl/dr_downloaders:"
    echo "    ./scripts/prepare_image.sh -i downloaders -d ccdl"
}

while getopts "phi:b:d:s:" opt; do
    case $opt in
        b)
            builder="--builder $OPTARG"
            ;;
        i)
            image=$OPTARG
            ;;
        d)
            dockerhub_repo=$OPTARG
            ;;
        p)
            pull="True"
            ;;
        s)
            service=$OPTARG
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

if [ -z "$image" ]; then
    echo "Error: you must specify an image with -i" >&2
    exit 1
fi

if [ -z "$service" ]; then
    service="workers"
fi

if [ -z "$dockerhub_repo" ]; then
    dockerhub_repo="ccdlstaging"
fi

# Default to "local" for system version if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="local$(date +%s)"
fi

# We want to check if a test image has been built for this branch. If
# it has we should use that rather than building it slowly.
IMAGE_NAME="$dockerhub_repo/dr_$image"
# shellcheck disable=SC2086
if [ "$(docker_img_exists $IMAGE_NAME $branch_name)" ] ; then
    docker pull "$IMAGE_NAME:$branch_name"
elif [ -n "$pull" ]; then
    docker pull "$IMAGE_NAME"
else
    echo ""
    echo "Rebuilding the $IMAGE_NAME image."
    finished=1
    attempts=0
    while [ $finished != 0 ] && [ $attempts -lt 1 ]; do
        if [ $attempts -gt 0 ]; then
            echo "Failed to build $IMAGE_NAME, trying again."
        fi

        CACHED_IMAGE="$IMAGE_NAME:latest"
        if test "$GITHUB_ACTIONS"; then
            # docker needs repositories to be lowercase
            CACHE_REPO="$(echo "ghrc.io/$GITHUB_REPOSITORY" | \
                tr '[:upper:]' '[:lower:]')"
            CACHED_IMAGE="$CACHE_REPO/dr_$image"
        fi

        DOCKER_BUILDKIT=1 docker buildx build $builder \
               -f "$service/dockerfiles/Dockerfile.$image" \
               -t "$IMAGE_NAME" \
               --build-arg BUILDKIT_INLINE_CACHE=1 \
               --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
               --cache-from $CACHED_IMAGE \
               --platform linux/amd64 \
               --progress plain \
               --push \
               .
        finished=$?
        attempts=$((attempts+1))
    done

    if [ $finished != 0 ] && [ $attempts -ge 3 ]; then
        echo "Could not build $IMAGE_NAME after three attempts."
        exit 1
    fi
fi
