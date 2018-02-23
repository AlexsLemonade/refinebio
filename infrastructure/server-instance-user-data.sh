#!/bin/bash

echo 'DOCKER_STORAGE_OPTIONS="--storage-driver overlay2"' > /etc/sysconfig/docker-storage
/etc/init.d/docker restart

git clone https://github.com/data-refinery/data-refinery.git

cd data-refinery && git fetch && git checkout origin terraform-docs

./install_nomad.sh

nomad agent -config infrastructure/server.hcl

# Currently doing dev just cause unsure how to handle prod secrets just yet.
./workers/format_nomad_with_env.sh

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 workers/downloader.nomad
nomad run -address http://$HOST_IP:4646 workers/processor.nomad
