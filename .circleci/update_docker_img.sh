# Build docker images
docker build -t ccdl/test_worker_base:$CIRCLE_TAG -f workers/Dockerfile.base .
# Build other images
# ...

# Log in and push images
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

# Push images
docker push ccdl/test_worker_base:$CIRCLE_TAG
# Push mnore images ...
