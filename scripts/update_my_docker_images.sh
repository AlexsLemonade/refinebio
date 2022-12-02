#!/bin/sh

# Exit on failure.
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of refinebio.
cd ..

print_description() {
    echo 'This script will re-build all refine.bio docker images and push '
    echo 'them to the specified Dockerhub repository.'
}

print_options() {
    cat << EOF
There are two required arguments for this script:
-d specifies the Dockerhub repository you would like to deploy to.
-v specifies the version you would like to build. This version will passed into
    the Docker image as the environment variable SYSTEM_VERSION.
    It also will be used as the tag for the Docker images built.

There is also optional arguments:
-a also build the affymetrix image
    (we normally don't because it is so intense to build)
-b sets a remote Docker builder for amd64 image builds.
-x specifies whether to build arm64 Docker images.
EOF
}

BUILDER=""

while getopts ":d:v:ab:xh" opt; do
    case $opt in
    d)
        export DOCKERHUB_REPO="$OPTARG"
        ;;
    v)
        export SYSTEM_VERSION="$OPTARG"
        ;;
    a)
        AFFYMETRIX=true
        ;;
    b)
        BUILDER="--builder $OPTARG"
        ;;
    x)  BUILD_ARM64=true
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
        echo "Option -$OPTARG requires an argument" >&2
        print_options >&2
        exit 1
        ;;
    esac
done

if [ -z "$DOCKERHUB_REPO" ]; then
    echo 'Error: must specify the Dockerhub repository with -d'
    exit 1
fi

if [ -z "$SYSTEM_VERSION" ]; then
    echo 'Error: must specify the version repository with -v'
    exit 1
fi

# Intentionally omit affymetrix unless specifically requested since it is so
# intense to build.
DOCKER_IMAGES="transcriptome smasher salmon no_op illumina downloaders compendia"
if [ "$AFFYMETRIX" ]; then
    DOCKER_IMAGES="$DOCKER_IMAGES affymetrix"
fi

# Set the version for the common project.
echo "$SYSTEM_VERSION" > common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist.
rm -f common/dist/*
(cd common && python3 setup.py sdist)  # Run in a subshell.

ARCH="$(uname -m)"

for DOCKER_IMAGE in $DOCKER_IMAGES; do
    case $DOCKER_IMAGE in
        api)
            DOCKER_FILE_PATH="api/dockerfiles/Dockerfile.api_production"
            ;;
        foreman)
            DOCKER_FILE_PATH="foreman/dockerfiles/Dockerfile.$DOCKER_IMAGE"
            ;;
        *)
            DOCKER_FILE_PATH="workers/dockerfiles/Dockerfile.$DOCKER_IMAGE"
            ;;
    esac

    IMAGE_NAME="$DOCKERHUB_REPO/dr_$DOCKER_IMAGE"

    echo "Building amd64 $IMAGE_NAME:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
    # shellcheck disable=SC2086
    DOCKER_BUILDKIT=1 docker buildx build $BUILDER \
        -f "$DOCKER_FILE_PATH" \
        -t "$IMAGE_NAME:$SYSTEM_VERSION" \
        -t "$IMAGE_NAME:latest" \
        --build-arg BUILDKIT_INLINE_CACHE=1 \
        --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
        --cache-from="$IMAGE_NAME:$SYSTEM_VERSION" \
        --cache-from="$IMAGE_NAME:latest" \
        --platform linux/amd64 \
        --push \
        .

    if [ "$ARCH" = "arm64" ] && [ "$BUILD_ARM64" ]; then
        echo "Building arm64 $IMAGE_NAME:$SYSTEM_VERSION image from $DOCKER_FILE_PATH."
        DOCKER_BUILDKIT=1 docker buildx build \
            -f "$DOCKER_FILE_PATH" \
            -t "$IMAGE_NAME:$SYSTEM_VERSION" \
            -t "$IMAGE_NAME:latest" \
            --build-arg BUILDKIT_INLINE_CACHE=1 \
            --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
            --cache-from="$IMAGE_NAME:$SYSTEM_VERSION" \
            --cache-from="$IMAGE_NAME:latest" \
            --platform linux/arm64 \
            --push \
            .
    fi
done
