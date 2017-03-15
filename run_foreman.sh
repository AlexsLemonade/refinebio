#!/bin/bash
HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run --link some-rabbit:rabbit --add-host=database:$HOST_IP test_master
