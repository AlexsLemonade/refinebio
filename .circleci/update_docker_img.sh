cd ~/refinebio

# Log into DockerHub
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

# Build docker images and push them into DockerHub
docker build -t ccdl/data_refinery_worker_base:$CIRCLE_TAG -f workers/Dockerfile.base .
docker push ccdl/data_refinery_worker_base:$CIRCLE_TAG

docker build -t ccdl/data_refinery_workers:$CIRCLE_TAG -f workers/Dockerfile .
docker push ccdl/data_refinery_workers:$CIRCLE_TAG

docker build -t ccdl/data_refinery_foreman:$CIRCLE_TAG -f foreman/Dockerfile .
docker push ccdl/data_refinery_foreman:$CIRCLE_TAG
