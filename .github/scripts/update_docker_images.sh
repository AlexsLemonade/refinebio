#!/bin/bash

set -e

# Load docker_image_exists function and $WORKER_IMAGES.
# shellcheck disable=SC1090
. ~/refinebio/scripts/common.sh

# Github won't set the branch name for us, so do it ourselves.
BRANCH=$(get_deploy_branch "$CI_TAG")
DOCKER_ACTION="--push"

if [[ "$BRANCH" == "master" ]]; then
    DOCKERHUB_REPO=ccdl
elif [[ "$BRANCH" == "dev" ]]; then
    DOCKERHUB_REPO=ccdlstaging
else
    echo "Why in the world was update_docker_images.sh called from a branch other than dev or master!?"
    exit 1
fi

echo "$CI_TAG" >~/refinebio/common/version

# Create ~/refinebio/common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist.
rm -f ~/refinebio/common/dist/*
cd ~/refinebio/common && python3 setup.py sdist

cd ~/refinebio

# Log into DockerHub.
docker login -u "$DOCKER_IO_USERNAME" -p "$DOCKER_IO_PASSWORD"

IMAGE_NAMES="$WORKER_IMAGES foreman api"
for IMAGE_NAME in $IMAGE_NAMES; do
    case $IMAGE_NAME in
    api)
        DOCKER_FILE_PATH="api/dockerfiles/Dockerfile.api_production"
        ;;
    foreman)
        DOCKER_FILE_PATH="foreman/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    *)
        DOCKER_FILE_PATH="workers/dockerfiles/Dockerfile.$IMAGE_NAME"
        ;;
    esac

    if docker_image_exists "$IMAGE_NAME" "$CI_TAG"; then
        echo "Docker image exists, skipping: $IMAGE_NAME:$CI_TAG"
    else
        echo "Building the $IMAGE_NAME:$CI_TAG image from $DOCKER_FILE_PATH."
        update_docker_image "$DOCKERHUB_REPO" "$IMAGE_NAME" "$CI_TAG" "$DOCKER_FILE_PATH" "$DOCKER_ACTION"
    fi
done
