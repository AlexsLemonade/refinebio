#!/bin/bash
HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run --link some-rabbit:rabbit --name worker1 --add-host=database:$HOST_IP -d bm_worker
