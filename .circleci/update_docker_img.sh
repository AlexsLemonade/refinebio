cd ~/refinebio
# Build docker images
docker build -t ccdl/data_refinery_worker_base:$CIRCLE_TAG -f workers/Dockerfile.base .
docker build -t ccdl/data_refinery_workers:$CIRCLE_TAG -f workers/Dockerfile .
docker build -t ccdl/data_refinery_foreman:$CIRCLE_TAG -f foreman/Dockerfile .

# Log into DockerHub and push images
docker login -u $DOCKER_ID -p $DOCKER_PASSWD
docker push ccdl/data_refinery_worker_base:$CIRCLE_TAG
docker push ccdl/data_refinery_workers:$CIRCLE_TAG
docker push ccdl/data_refinery_foreman:$CIRCLE_TAG
