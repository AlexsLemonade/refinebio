#!/bin/sh

set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# shellcheck disable=SC1091
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

DOCKER_ACTION="--push"

# Intentionally omit affymetrix unless specifically requested since it is so
# intense to build.
IMAGE_NAMES="base migrations common_tests foreman api_base api api_local \
    transcriptome smasher salmon no_op illumina downloaders compendia"
if [ "$BUILD_AFFYMETRIX" ]; then
    IMAGE_NAMES="$IMAGE_NAMES affymetrix"
fi

# Set the version for the common project.
echo "$SYSTEM_VERSION" >common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
# Remove old common distributions if they exist.
rm -f common/dist/*
(cd common && python3 setup.py sdist 1>/dev/null) # Run quietly in a subshell.

for IMAGE_NAME in $IMAGE_NAMES; do
    case $IMAGE_NAME in
    api_base | api | api_local)
        DOCKER_FILE_PATH="api/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    base | common_tests | migrations)
        DOCKER_FILE_PATH="common/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    foreman)
        DOCKER_FILE_PATH="foreman/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    *)
        DOCKER_FILE_PATH="workers/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    esac

    echo "Building the $IMAGE_NAME:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
    update_docker_image "$DOCKERHUB_REPO" "$IMAGE_NAME" "$SYSTEM_VERSION" "$DOCKER_FILE_PATH" "$DOCKER_ACTION"
done
