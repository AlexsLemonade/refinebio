#!/bin/bash

if [[ -z $DOCKER_ID ]]; then
    echo 'Please set $DOCKER_ID!'
    exit 1
fi

# Log into DockerHub
echo "Please enter your password for Dockerhub.com"
docker login -u $DOCKER_ID

if [[ -z $DOCKERHUB_REPO ]]; then
    DOCKERHUB_REPO=$DOCKER_ID
fi

# Intentionally omit affymetrix since it is so intense to build.
CCDL_WORKER_IMGS="salmon transcriptome no_op downloaders illumina"

for IMG in $CCDL_WORKER_IMGS; do
    image_name=$DOCKERHUB_REPO/dr_$IMG
    # Build and push image
    docker build -t $image_name -f workers/dockerfiles/Dockerfile.$IMG .
    docker push $image_name
done

# Build and push foreman image
docker build -t $DOCKERHUB_REPO/dr_foreman -f foreman/dockerfiles/Dockerfile.foreman .
docker push $DOCKERHUB_REPO/dr_foreman

# Build and push API image
docker build -t $DOCKERHUB_REPO/dr_api -f api/dockerfiles/Dockerfile.api_production .
docker push $DOCKERHUB_REPO/dr_api
