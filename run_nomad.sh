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
      > nomad.logs &

docker build -t dr_worker -f workers/Dockerfile .

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

if [ ! -f workers/downloader.nomad -o ! -f workers/processor.nomad ]; then
    echo "One or more Job Specifications for Nomad are missing, building them from the template."
    ./workers/format_nomad_with_env.sh
fi

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 workers/downloader.nomad
nomad run -address http://$HOST_IP:4646 workers/processor.nomad
