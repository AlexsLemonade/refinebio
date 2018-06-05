#!/bin/bash

if [[ -z $DOCKER_ID ]]; then
    echo 'Please set $DOCKER_ID!'
    exit 1
fi

# Log into DockerHub
echo "Please enter your password for Dockerhub.com"
docker login -u $DOCKER_ID

# Intentionally omit affymetrix since it is so intense to build.
CCDL_WORKER_IMGS="salmon transcriptome no_op downloaders illumina"

cd ~/refinebio
for IMG in $CCDL_WORKER_IMGS; do
    image_name=$DOCKER_ID/dr_$IMG
    # Build and push image
    echo "Building $image_name"
    docker build -t $image_name -f workers/dockerfiles/Dockerfile.$IMG .
    docker push $image_name
done

# Build and push foreman image
echo "Building $image_name"
docker build -t $DOCKER_ID/data_refinery_foreman -f foreman/dockerfiles/Dockerfile.foreman .
docker push $DOCKER_ID/data_refinery_foreman

# Build and push API image
echo "Building $image_name"
docker build -t $DOCKER_ID/data_refinery_api -f api/dockerfiles/Dockerfile.api_production .
docker push $DOCKER_ID/data_refinery_api
