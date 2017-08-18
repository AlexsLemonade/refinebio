#!/bin/bash
echo ECS_CLUSTER=data-refinery >> /etc/ecs/ecs.config

echo 'DOCKER_STORAGE_OPTIONS="--storage-driver overlay2"' > /etc/sysconfig/docker-storage
/etc/init.d/docker restart
