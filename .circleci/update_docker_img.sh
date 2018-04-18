# This function checks whether a given docker image name ($1:$CIRCLE_TAG)
# exists in Docker Hub or not using Docker Hub API V2. Based on:
# https://stackoverflow.com/questions/32113330/check-if-imagetag-combination-already-exists-on-docker-hub
function docker_img_exists() {
    TOKEN=$(curl -s -H "Content-Type: application/json" -X POST \
                 -d '{"username": "'${DOCKER_ID}'", "password": "'${DOCKER_PASSWD}'"}' \
                 https://hub.docker.com/v2/users/login/ | jq -r .token)
    EXISTS=$(curl -s -H "Authorization: JWT ${TOKEN}" \
                  https://hub.docker.com/v2/repositories/$1/tags/?page_size=10000 \
             | jq -r "[.results | .[] | .name == \"${CIRCLE_TAG}\"] | any")
    test $EXISTS = true
}

# Docker images that we want to protect from accidental overwriting in "ccdl" account.
CCDL_IMGS="data_refinery_worker_base data_refinery_workers data_refinery_foreman"

# If any of the three images could be overwritten by the building process,
# print out a message and terminate immediately.
for IMG in $CCDL_IMGS; do
    if docker_img_exists ccdl/$IMG; then
        echo "Docker image exists, building process terminated: $IMG:$CIRCLE_TAG"
        exit
    fi
done

# Create ~/refinebio/common/dist/data-refinery-common-*.tar.gz, which is
# required by data_refinery_workers and data_refinery_foreman images.
cd ~/refinebio/common && python setup.py sdist

# Log into DockerHub
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

cd ~/refinebio
# Build and push worker_base image
docker build -t ccdl/data_refinery_worker_base:$CIRCLE_TAG -f workers/Dockerfile.base .
docker push ccdl/data_refinery_worker_base:$CIRCLE_TAG
# Update latest version
docker tag ccdl/data_refinery_worker_base:$CIRCLE_TAG ccdl/data_refinery_worker_base:latest
docker push ccdl/data_refinery_worker_base:latest

# Build and push workers image
docker build -t ccdl/data_refinery_workers:$CIRCLE_TAG -f workers/Dockerfile .
docker push ccdl/data_refinery_workers:$CIRCLE_TAG
# Update latest version
docker tag ccdl/data_refinery_workers:$CIRCLE_TAG ccdl/data_refinery_workers:latest
docker push ccdl/data_refinery_workers:latest

# Build and push foreman image
docker build -t ccdl/data_refinery_foreman:$CIRCLE_TAG -f foreman/Dockerfile .
docker push ccdl/data_refinery_foreman:$CIRCLE_TAG
# Update latest version
docker tag ccdl/data_refinery_foreman:$CIRCLE_TAG ccdl/data_refinery_foreman:latest
docker push ccdl/data_refinery_foreman:latest
