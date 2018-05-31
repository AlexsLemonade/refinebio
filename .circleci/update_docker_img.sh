#!/bin/bash

# Load docker_img_exists function
source ~/refinebio/common.sh

# Docker images that we want to protect from accidental overwriting in "ccdl" account.
CCDL_WORKER_IMGS="illumina affymetrix salmon transcriptome no_op downloaders"

# If any of the three images could be overwritten by the building process,
# print out a message and terminate immediately.
for IMG in $CCDL_WORKER_IMGS; do
    image_name=ccdl/dr_$IMG
    if docker_img_exists $image_name $CIRCLE_TAG; then
        echo "Docker image exists, building process terminated: $image_name:$CIRCLE_TAG"
        exit
    fi
done

CCDL_OTHER_IMGS="ccdl/data_refinery_foreman ccdl/data_refinery_api"

# Handle the foreman separately.
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
for IMG in $CCDL_WORKER_IMGS; do
    image_name=ccdl/dr_$IMG
    # Build and push image
    docker build -t $image_name:$CIRCLE_TAG -f workers/dockerfiles/$IMG .
    docker push $image_name:$CIRCLE_TAG
    # Update latest version
    docker tag $image_name:$CIRCLE_TAG $image_name:latest
    docker push $image_name:latest
done

# Build and push foreman image
./prepare_image.sh -i foreman -s foreman
docker push ccdl/data_refinery_foreman:$CIRCLE_TAG
# Update latest version
docker tag ccdl/data_refinery_foreman:$CIRCLE_TAG ccdl/data_refinery_foreman:latest
docker push ccdl/data_refinery_foreman:latest

# Build and push API image
./prepare_image.sh -i api_production -s api
docker tag ccdl/data_refinery_api ccdl/data_refinery_api:$CIRCLE_TAG
docker push ccdl/data_refinery_api:$CIRCLE_TAG
# Update latest version
docker tag ccdl/data_refinery_api:$CIRCLE_TAG ccdl/data_refinery_api:latest
docker push ccdl/data_refinery_api:latest
