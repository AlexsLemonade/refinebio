#!/bin/bash

# This script runs the Data Refinery Workers locally.
# What that actually means is that it runs Nomad and registers the
# Nomad jobs for the Data Refinery Workers. Nomad is listening for
# messages telling it to run Processor and Downloader jobs. When it
# receives them, it will run the jobs within new Docker containers.

# Figure out the right location to put the nomad directory.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

source common.sh
HOST_IP=$(get_ip_address)

echo "Rebuilding Docker image while waiting for Nomad to come online."

nomad_dir="$script_directory/nomad_dir"
if [ ! -d $nomad_dir ]; then
    mkdir nomad_dir
fi

# Start the nomad in development locally
nomad agent -bind $HOST_IP \
      -data-dir $nomad_dir \
      -dev \
      &> nomad.logs &

# While we wait for Nomad to start, make sure the Docker registry has
# an up-to-date Docker image for the workers sub-project. We run a
# local Docker registry because Nomad must pull from a registry, it
# will not use an image which lives on the host machine.

# Make sure the local Docker registry is running.
if [[ -z $(docker ps | grep "image-registry") ]]; then
    docker run -d -p 5000:5000 --restart=always --name image-registry registry:2
fi

# Build, tag, and push an image for the workers to the local registry.
docker build -t dr_worker -f workers/Dockerfile .
docker tag dr_worker localhost:5000/dr_worker
docker push localhost:5000/dr_worker

# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  http://$HOST_IP:4646/v1/status/leader)
}

# Wait for Nomad to get started.
nomad_status=$(check_nomad_status)
while [ $nomad_status != "200" ]; do
    sleep 1
    nomad_status=$(check_nomad_status)
done

echo "Nomad is online. Registering jobs."

./workers/format_nomad_with_env.sh

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 workers/downloader.nomad
nomad run -address http://$HOST_IP:4646 workers/processor.nomad
