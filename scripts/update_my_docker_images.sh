#!/bin/sh

# This script is very similar to .circleci/update_docker_images.sh but it has less
# production/cloud related checks.

# Exit on failure
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of refinebio
cd ..

print_description() {
    echo 'This script will re-build all refine.bio docker images and push them to'
    echo 'the specified Dockerhub repo.'
}

print_options() {
    cat << EOF
There are two required arguments for this script:
-d specifies the Dockerhub repo you would like to deploy to.
-v specifies the version you would like to build. This version will passed into
    the Docker image as the environment variable SYSTEM_VERSION.
    It also will be used as the tag for the Docker images built.

There is also one optional argument:
-a also build the affymetrix image
    (we normally don't because it is so intense to build)
EOF
}

while getopts ":d:v:ah" opt; do
    case $opt in
    d)
        export DOCKERHUB_REPO=$OPTARG
        ;;
    v)
        export SYSTEM_VERSION=$OPTARG
        ;;
    a)
        AFFYMETRIX=true
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
    echo 'Error: must specify the Dockerhub repo with -d'
    exit 1
fi

if [ -z "$SYSTEM_VERSION" ]; then
    echo 'Error: must specify the version repo with -v'
    exit 1
fi

# Intentionally omit affymetrix unless specifically requested since it is so intense to build.
CCDL_WORKER_IMGS="salmon transcriptome no_op downloaders illumina smasher"
if [ "$AFFYMETRIX" ]; then
    CCDL_WORKER_IMGS="$CCDL_WORKER_IMGS affymetrix"
fi

# Set the version for the common project.
echo "$SYSTEM_VERSION" > common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist
rm -f common/dist/*
(cd common && python3 setup.py sdist)

for IMG in $CCDL_WORKER_IMGS; do
    image_name="$DOCKERHUB_REPO/dr_$IMG"

    echo "Building docker image: $image_name:$SYSTEM_VERSION"
    # Build and push image.
    docker build \
           -t "$image_name:$SYSTEM_VERSION" \
           -f "workers/dockerfiles/Dockerfile.$IMG" \
           --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" .
    docker push "$image_name:$SYSTEM_VERSION"
    # Update latest version
    docker tag "$image_name:$SYSTEM_VERSION" "$image_name:latest"
    docker push "$image_name:latest"
done

# Build and push Foreman image.
FOREMAN_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_foreman"
docker build \
       -t "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION" \
       -f foreman/dockerfiles/Dockerfile.foreman \
       --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" .
docker push "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION"
# Update latest version
docker tag "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION" "$FOREMAN_DOCKER_IMAGE:latest"
docker push "$FOREMAN_DOCKER_IMAGE:latest"

# Build and push API image.
API_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_api"
docker build \
       -t "$API_DOCKER_IMAGE:$SYSTEM_VERSION" \
       -f api/dockerfiles/Dockerfile.api_production \
       --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" .
docker push "$API_DOCKER_IMAGE:$SYSTEM_VERSION"
# Update latest version
docker tag "$API_DOCKER_IMAGE:$SYSTEM_VERSION" "$API_DOCKER_IMAGE:latest"
docker push "$API_DOCKER_IMAGE:latest"
