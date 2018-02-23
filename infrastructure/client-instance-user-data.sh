#!/bin/bash

echo 'DOCKER_STORAGE_OPTIONS="--storage-driver overlay2"' > /etc/sysconfig/docker-storage
/etc/init.d/docker restart

git clone https://github.com/data-refinery/data-refinery.git

cd data-refinery && git fetch && git checkout origin terraform-docs

./install_nomad.sh

nomad agent -config infrastructure/client.hcl
