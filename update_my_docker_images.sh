#!/bin/bash -e

# This script is very similar to .circleci/update_docker_images.sh but it has less
# production/cloud related checks.

print_description() {
    echo 'This script will re-build all refine.bio docker images and push them to'
    echo 'the specified Dockerhub repo.'
}

print_options() {
    echo 'There is are two arguments for this script: -d and -v, and they are not optional.'
    echo '-d specifies the Dockerhub repo you would like to deploy to.'
    echo '-v specifies the version you would like to build. This version will passed into'
    echo '  the Docker image as the environment variable SYSTEM_VERSION.'
    echo '  It also will be used as the tag for the Docker images built.'
}

while getopts ":d:v:h" opt; do
    case $opt in
    d)
        export DOCKERHUB_REPO=$OPTARG
        ;;
    v)
        export SYSTEM_VERSION=$OPTARG
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

if [[ -z $DOCKERHUB_REPO ]]; then
    echo 'Error: must specify the Dockerhub repo with -d'
    exit 1
fi

if [[ -z $SYSTEM_VERSION ]]; then
    echo 'Error: must specify the version repo with -v'
    exit 1
fi

# Intentionally omit affymetrix since it is so intense to build.
CCDL_WORKER_IMGS="salmon transcriptome no_op downloaders illumina smasher"

# Set the version for the common project.
echo $SYSTEM_VERSION > common/version

# Create common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist
rm -f common/dist/*
cd common && python setup.py sdist

cd ..
for IMG in $CCDL_WORKER_IMGS; do
    image_name="$DOCKERHUB_REPO/dr_$IMG"

    echo "Building docker image: $image_name:$SYSTEM_VERSION"
    # Build and push image.
    docker build \
           -t $image_name:$SYSTEM_VERSION \
           -f workers/dockerfiles/Dockerfile.$IMG \
           --build-arg SYSTEM_VERSION=$SYSTEM_VERSION .
    docker push $image_name:$SYSTEM_VERSION
    # Update latest version
    docker tag $image_name:$SYSTEM_VERSION $image_name:latest
    docker push $image_name:latest
done

# Build and push Foreman image.
FOREMAN_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_foreman"
docker build \
       -t "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION" \
       -f foreman/dockerfiles/Dockerfile.foreman \
       --build-arg SYSTEM_VERSION=$SYSTEM_VERSION .
docker push "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION"
# Update latest version
docker tag "$FOREMAN_DOCKER_IMAGE:$SYSTEM_VERSION" "$FOREMAN_DOCKER_IMAGE:latest"
docker push "$FOREMAN_DOCKER_IMAGE:latest"

# Build and push API image.
API_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_api"
docker build \
       -t "$API_DOCKER_IMAGE:$SYSTEM_VERSION" \
       -f api/dockerfiles/Dockerfile.api_production \
       --build-arg SYSTEM_VERSION=$SYSTEM_VERSION .
docker push "$API_DOCKER_IMAGE:$SYSTEM_VERSION"
# Update latest version
docker tag "$API_DOCKER_IMAGE:$SYSTEM_VERSION" "$API_DOCKER_IMAGE:latest"
docker push "$API_DOCKER_IMAGE:latest"
