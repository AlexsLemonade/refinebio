#!/bin/sh

# Exit on failure.
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

. ./common.sh

# Get access to all of refinebio.
cd ..

print_description() {
    echo 'This script will re-build all refine.bio docker images and push '
    echo 'them to the specified Dockerhub repository.'
}

print_options() {
    cat <<EOF
There are two required arguments for this script:
-r specifies the Dockerhub repository you would like to deploy to.
-v specifies the version you would like to build. This version will passed into
    the Docker image as the environment variable SYSTEM_VERSION.
    It also will be used as the tag for the Docker images built.

There are also optional arguments:
-a also build the affymetrix image
    (we normally don't because it is so intense to build)
EOF
}

while getopts ":r:v:ah" opt; do
    case $opt in
    r)
        export DOCKERHUB_REPO="$OPTARG"
        ;;
    v)
        export SYSTEM_VERSION="$OPTARG"
        ;;
    a)
        BUILD_AFFYMETRIX=true
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

if [ -z "$DOCKERHUB_REPO" ]; then
    echo 'Error: must specify the Dockerhub repository with -r'
    exit 1
fi

if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="$(get_branch_hash)"
fi

# Intentionally omit affymetrix unless specifically requested since it is so
# intense to build.
image_names="base migrations common_tests foreman api_base api_production api_local \
    transcriptome smasher salmon no_op illumina downloaders compendia"
if [ "$BUILD_AFFYMETRIX" ]; then
    image_names="$image_names affymetrix"
fi

# Set the version for the common project.
echo "$SYSTEM_VERSION" >common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
# Remove old common distributions if they exist.
rm -f common/dist/*
(cd common && python3 setup.py sdist 1>/dev/null) # Run quietly in a subshell.

# shellcheck disable=SC2086
for image_name in $image_names; do
    case $image_name in
    api_base | api_local | api_production)
        DOCKER_FILE_PATH="api/dockerfiles/Dockerfile.$image_name"
        ;;
    base | common_tests | migrations)
        DOCKER_FILE_PATH="common/dockerfiles/Dockerfile.$image_name"
        ;;
    foreman)
        DOCKER_FILE_PATH="foreman/dockerfiles/Dockerfile.$image_name"
        ;;
    *)
        DOCKER_FILE_PATH="workers/dockerfiles/Dockerfile.$image_name"
        ;;
    esac

    echo
    echo "Building the $image_name:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
    echo

    DOCKERHUB_IMAGE="$DOCKERHUB_REPO/dr_$image_name"
    CACHE_FROM_LATEST="cache-from=type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest"
    CACHE_FROM_VERSION="cache-from=type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION"
    CACHE_TO_LATEST="cache-to=type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest,mode=max"
    CACHE_TO_VERSION="cache-to=type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION,mode=max"

    set_up_docker_builder

    docker buildx build \
        --build-arg DOCKERHUB_REPO="$DOCKERHUB_REPO" \
        --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
        --"$CACHE_FROM_LATEST" \
        --"$CACHE_FROM_VERSION" \
        --"$CACHE_TO_LATEST" \
        --"$CACHE_TO_VERSION" \
        --file "$DOCKER_FILE_PATH" \
        --platform linux/amd64 \
        --push \
        --tag "$DOCKERHUB_IMAGE:latest" \
        --tag "$DOCKERHUB_IMAGE:$SYSTEM_VERSION" \
        .
done
