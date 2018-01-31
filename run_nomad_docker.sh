#! /bin/bash

# This script starts a Nomad agent running within a container
# locally. In order to do this and have it work, we need to bind the
# nomad agent to an IP address other than localhost. To do this, we
# need to be able to specify the IP address of the container. Docker
# only lets you specify the IP address of a container if you do so
# within a subnet you have created. Therefore we first create a
# subnet, and specify that subnet for all containers which need to
# interact with it.

# If data-refinery Docker network doesn't exist, create it.
if [[ -z $(docker network ls | grep data-refinery) ]]; then
    # The subnet 172.29.0.0/16 was chosen because it is within
    # Docker's range of 172.0.0.0/24 and doesn't seem like it will
    # have conflicts with other Docker services.
    docker network create --subnet=172.29.0.0/16 data-refinery
fi

# If there is a stopped "nomad" container, remove it so the dev
# doesn't have to prune it herself.
docker container rm nomad

# 172.29.0.1 is reserved for Docker itself, so use the next IP after
# that.
container_ip=172.29.0.2

docker run -d \
  --net data-refinery \
  --ip $container_ip \
  --name nomad \
  -p 4646:4646 \
  -p 4647:4647 \
  -p 4648:4648 \
  miserlou/nomad agent -dev -bind $container_ip -network-interface eth0
