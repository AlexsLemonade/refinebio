#!/bin/bash
echo ECS_CLUSTER=data-refinery >> /etc/ecs/ecs.config

echo 'OPTIONS="${OPTIONS} --storage-opt dm.basesize=40G"' >> /etc/sysconfig/docker
/etc/init.d/docker restart
