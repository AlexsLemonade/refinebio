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
-b sets a remote Docker builder.
-x specifies target Docker images architecture (amd64|arm64).
EOF
}

while getopts ":r:v:abx:h" opt; do
    case $opt in
    r)
        export DOCKERHUB_REPO="$OPTARG"
        ;;
    v)
        export SYSTEM_VERSION="$OPTARG"
        ;;
    a)
        AFFYMETRIX=true
        ;;
    b)
        DOCKER_BUILDER="$OPTARG"
        ;;
    x)
        SYSTEM_ARCH="$OPTARG"
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

if [ -z "$SYSTEM_ARCH" ]; then
    SYSTEM_ARCH="amd64"
fi

if [ "$SYSTEM_ARCH" != "amd64" ] && [ "$SYSTEM_ARCH" != "arm64" ]; then
    echo "Error: unsupported architecture $SYSTEM_ARCH: use either amd64 or arm64."
    exit 1
fi

if [ -z "$DOCKERHUB_REPO" ]; then
    echo 'Error: must specify the Dockerhub repository with -r'
    exit 1
fi

if [ -z "$SYSTEM_VERSION" ]; then
    echo 'Error: must specify the version repository with -v'
    exit 1
fi

if [ -z "$DOCKER_BUILDER" ]; then
    DOCKER_BUILDER="local"
fi


# Intentionally omit affymetrix unless specifically requested since it is so
# intense to build.
DOCKER_IMAGES="transcriptome smasher salmon no_op illumina downloaders compendia"
if [ "$AFFYMETRIX" ]; then
    DOCKER_IMAGES="$DOCKER_IMAGES affymetrix"
fi

# Set the version for the common project.
echo "$SYSTEM_VERSION" >common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist.
rm -f common/dist/*
(cd common && python3 setup.py sdist 1>/dev/null) # Run quietly in a subshell.

ARCH="$(uname -m)"

if [ "$ARCH" != "$SYSTEM_ARCH" ] && [ -z "$DOCKER_BUILDER" ]; then
    echo "You're building $SYSTEM_ARCH images on the $ARCH machine without " \
        "specifying a DOCKER_BUILDER. This may take a while..."
    sleep 5
fi

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

    echo "Building $SYSTEM_ARCH $IMAGE_NAME:$SYSTEM_VERSION image from " \
        "$DOCKER_FILE_PATH."

    if test "$DOCKER_BUILDER"; then
        echo "Using builder $DOCKER_BUILDER."
        docker buildx use "$DOCKER_BUILDER"
    fi

    DOCKER_BUILDKIT=1 docker buildx build \
        --build-arg BUILDKIT_INLINE_CACHE=1 \
        --build-arg DOCKERHUB_REPO="$DOCKERHUB_REPO" \
        --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
        --cache-from "$IMAGE_NAME:$SYSTEM_VERSION" \
        --cache-from "$IMAGE_NAME:latest" \
        --file "$DOCKER_FILE_PATH" \
        --platform "linux/$SYSTEM_ARCH" \
        --push \
        --tag "$IMAGE_NAME:$SYSTEM_VERSION" \
        --tag "$IMAGE_NAME:latest" \
        .

    docker pull \
        --platform "linux/$SYSTEM_ARCH" \
        "$IMAGE_NAME:latest"
done
