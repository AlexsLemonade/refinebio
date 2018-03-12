CIRCLE_TAG=dhutest # dhu test only

cd ~/refinebio
# Build docker images
docker build -t ccdl/test_worker_base:$CIRCLE_TAG -f workers/Dockerfile.base .
docker build -t ccdl/test_workers:$CIRCLE_TAG -f workers/Dockerfile .
docker build -t ccdl/test_foreman:$CIRCLE_TAG -f foreman/Dockerfile .


# Log in and push images
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

# Push images
docker push ccdl/test_worker_base:$CIRCLE_TAG
docker push ccdl/test_worker:$CIRCLE_TAG
docker push ccdl/test_foreman:$CIRCLE_TAG
