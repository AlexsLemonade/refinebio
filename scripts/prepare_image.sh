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
    echo
    echo "Examples:"
    echo "    Build the image ccdl/dr_downloaders:"
    echo "    ./scripts/prepare_image.sh -i downloaders -d ccdl"
}

while getopts "phi:d:s:" opt; do
    case $opt in
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
image_name="$dockerhub_repo/dr_$image"
# shellcheck disable=SC2086
if [ "$(docker_img_exists $image_name $branch_name)" ] ; then
    docker pull "$image_name:$branch_name"
elif [ -n "$pull" ]; then
    docker pull "$image_name"
else
    echo ""
    echo "Rebuilding the $image_name image."
    finished=1
    attempts=0
    while [ $finished != 0 ] && [ $attempts -lt 3 ]; do
        if [ $attempts -gt 0 ]; then
            echo "Failed to build $image_name, trying again."
        fi


        if $GITHUB_ACTIONS; then
            CACHE_REPO="docker.pkg.github.com/$GITHUB_REPOSITORY"
            CACHED_PACKAGE="$CACHE_REPO/dr_$image"
            CACHE="--build-arg BUILDKIT_INLINE_CACHE=1 --cache-from $CACHED_PACKAGE"
        fi

        docker build \
               -t "$image_name" \
               -f "$service/dockerfiles/Dockerfile.$image" \
               --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
               $CACHE .
        finished=$?
        attempts=$((attempts+1))
    done

    if [ $finished != 0 ] && [ $attempts -ge 3 ]; then
        echo "Could not build $image_name after three attempts."
        exit 1
    fi
fi
