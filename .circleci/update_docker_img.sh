# Build docker images
docker build -t dongbohu/worker_base -f workers/Dockerfile.base .
# Build other images
# ...

# Log in and push images
#echo $docker_passwd | docker login -u $docker_id --password-std
docker login -u $docker_id -p $docker_passwd

# Push images
docker push dongbohu/worker_base
# Push mnore images ...
