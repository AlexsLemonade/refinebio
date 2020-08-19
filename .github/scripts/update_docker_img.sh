#!/bin/bash
set -e

# Load docker_img_exists function and $CCDL_WORKER_IMAGES
source ~/refinebio/scripts/common.sh

# Github won't set the branch name for us, so do it ourselves.
branch=$(get_master_or_dev "$CI_TAG")

if [[ "$branch" == "master" ]]; then
    DOCKERHUB_REPO=ccdl
elif [[ "$branch" == "dev" ]]; then
    DOCKERHUB_REPO=ccdlstaging
else
    echo "Why in the world was update_docker_img.sh called from a branch other than dev or master?!?!?"
    exit 1
fi

echo "$CI_TAG" > ~/refinebio/common/version

# Create ~/refinebio/common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
## Remove old common distributions if they exist
rm -f ~/refinebio/common/dist/*
cd ~/refinebio/common && python3 setup.py sdist

# Log into DockerHub
docker login -u "$DOCKER_ID" -p "$DOCKER_PASSWD"

cd ~/refinebio
for IMAGE in $CCDL_WORKER_IMAGES; do
    image_name="$DOCKERHUB_REPO/dr_$IMAGE"
    if docker_img_exists "$image_name" "$CI_TAG"; then
        echo "Docker image exists, skipping: $image_name:$CI_TAG"
    else
        echo "Building docker image: $image_name:$CI_TAG"
        # Build and push image. We use the CI_TAG as the system version.
        docker build \
               -t "$image_name:$CI_TAG" \
               -f "workers/dockerfiles/Dockerfile.$IMAGE" \
               --build-arg SYSTEM_VERSION="$CI_TAG" .
        docker push "$image_name:$CI_TAG"
        # Update latest version
        docker tag "$image_name:$CI_TAG" "$image_name:latest"
        docker push "$image_name:latest"

        # Save some space when we're through
        docker rmi "$image_name:$CI_TAG"
    fi
done

# Build and push foreman image
FOREMAN_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_foreman"
if docker_img_exists "$FOREMAN_DOCKER_IMAGE" "$CI_TAG"; then
    echo "Docker image exists, skipping: $FOREMAN_DOCKER_IMAGE:$CI_TAG"
else
    # Build and push image. We use the CI_TAG as the system version.
    docker build \
           -t "$FOREMAN_DOCKER_IMAGE:$CI_TAG" \
           -f foreman/dockerfiles/Dockerfile.foreman \
           --build-arg SYSTEM_VERSION="$CI_TAG" .
    docker push "$FOREMAN_DOCKER_IMAGE:$CI_TAG"
    # Update latest version
    docker tag "$FOREMAN_DOCKER_IMAGE:$CI_TAG" "$FOREMAN_DOCKER_IMAGE:latest"
    docker push "$FOREMAN_DOCKER_IMAGE:latest"
fi

# Build and push API image
API_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_api"
if docker_img_exists "$API_DOCKER_IMAGE" "$CI_TAG"; then
    echo "Docker image exists, skipping: $API_DOCKER_IMAGE:$CI_TAG"
else
    # Build and push image. We use the CI_TAG as the system version.
    docker build \
           -t "$API_DOCKER_IMAGE:$CI_TAG" \
           -f api/dockerfiles/Dockerfile.api_production \
           --build-arg SYSTEM_VERSION="$CI_TAG" .
    docker push "$API_DOCKER_IMAGE:$CI_TAG"
    # Update latest version
    docker tag "$API_DOCKER_IMAGE:$CI_TAG" "$API_DOCKER_IMAGE:latest"
    docker push "$API_DOCKER_IMAGE:latest"
fi
