#! /bin/bash

# If data-refinery Docker network doesn't exist, create it.
if [[ -z $(docker network ls | grep data-refinery) ]]; then
    docker network create --subnet=172.29.0.0/16 data-refinery
fi

docker run -d \
  --net data-refinery \
  --ip 172.29.0.2 \
  --name nomad \
  -p 4646:4646 \
  -p 4647:4647 \
  -p 4648:4648 \
  miserlou/nomad agent -dev -bind 172.29.0.2 -network-interface eth0
