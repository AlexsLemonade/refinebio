#!/bin/bash

# This won't be testable without this var, so mock it out
CIRCLE_TAG="vTesting"

# Load docker_img_exists function
source ~/refinebio/common.sh

if [ $CIRCLE_BRANCH == "master" ]; then
    DOCKERHUB_REPO=ccdl
elif [[ $CIRCLE_BRANCH == "dev" || $CIRCLE_BRANCH == "kurtwheeler-addtional-envs" ]]; then
    DOCKERHUB_REPO=ccdlstaging
else
    echo "Why in the world was update_docker_img.sh called from a branch other than `dev` or `master`?!?!?"
    exit 1
fi

# Docker images that we want to protect from accidental overwriting in "ccdl" account.
CCDL_WORKER_IMGS="illumina affymetrix salmon transcriptome no_op downloaders"

# If any of the three images could be overwritten by the building process,
# print out a message and terminate immediately.
for IMG in $CCDL_WORKER_IMGS; do
    image_name="$DOCKERHUB_REPO/dr_$IMG"
    if docker_img_exists $image_name $CIRCLE_TAG; then
        echo "Docker image exists, building process terminated: $image_name:$CIRCLE_TAG"
        exit
    fi
done

CCDL_OTHER_IMGS="$DOCKERHUB_REPO/dr_foreman $DOCKERHUB_REPO/dr_api"

# Handle the foreman and API separately.
for IMG in $CCDL_OTHER_IMGS; do
    if docker_img_exists $IMG $CIRCLE_TAG; then
        echo "Docker image exists, building process terminated: $IMG:$CIRCLE_TAG"
        exit
    fi
done

# Create ~/refinebio/common/dist/data-refinery-common-*.tar.gz, which is
# required by the workers and data_refinery_foreman images.
cd ~/refinebio/common && python setup.py sdist

# Log into DockerHub
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

cd ~/refinebio
# for IMG in $CCDL_WORKER_IMGS; do
#     image_name="$DOCKERHUB_REPO/dr_$IMG"
#     # Build and push image
#     docker build -t $image_name:$CIRCLE_TAG -f workers/dockerfiles/Dockerfile.$IMG .
#     docker push $image_name:$CIRCLE_TAG
#     # Update latest version
#     docker tag $image_name:$CIRCLE_TAG $image_name:latest
#     docker push $image_name:latest
# done

# Build and push foreman image
FOREMAN_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_foreman"
docker build -t "$FOREMAN_DOCKER_IMAGE:$CIRCLE_TAG" -f foreman/dockerfiles.foreman .
docker push "$FOREMAN_DOCKER_IMAGE:$CIRCLE_TAG"
# Update latest version
docker tag "$FOREMAN_DOCKER_IMAGE:$CIRCLE_TAG" "$FOREMAN_DOCKER_IMAGE:latest"
docker push "$FOREMAN_DOCKER_IMAGE:latest"

# Build and push API image
API_DOCKER_IMAGE="$DOCKERHUB_REPO/dr_api"
# docker tag "$API_DOCKER_IMAGE" "$API_DOCKER_IMAGE:$CIRCLE_TAG"
docker build -t "$API_DOCKER_IMAGE:$CIRCLE_TAG" -f api/dockerfiles/Dockerfile.api_production .
docker push "$API_DOCKER_IMAGE:$CIRCLE_TAG"
# Update latest version
docker tag "$API_DOCKER_IMAGE:$CIRCLE_TAG" "$API_DOCKER_IMAGE:latest"
docker push "$API_DOCKER_IMAGE:latest"
