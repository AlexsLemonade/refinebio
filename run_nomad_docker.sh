#! /bin/bash
docker run -d \
  --name nomad \
  -p 4646:4646 \
  -p 4647:4647 \
  -p 4648:4648 \
  miserlou/nomad agent -dev -network-interface eth0
