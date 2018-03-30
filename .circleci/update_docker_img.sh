# Create ~/refinebio/common/dist/data-refinery-common-*.tar.gz, which is
# required by data_refinery_workers and data_refinery_foreman images.
cd ~/refinebio/common && python setup.py sdist

# Log into DockerHub
cd ~/refinebio
docker login -u $DOCKER_ID -p $DOCKER_PASSWD

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
