#!/bin/bash

# This script runs the Data Refinery Workers locally.
# What that actually means is that it runs Nomad and registers the
# Nomad jobs for the Data Refinery Workers. Nomad is listening for
# messages telling it to run Processor and Downloader jobs. When it
# receives them, it will run the jobs within new Docker containers.

print_description() {
    echo "Starts Nomad and registers the jobs with it. This involves re-building the "
    echo "Docker images and running format_nomad_with_env.sh to format the Nomad job "
    echo "specifications for the correct environment."
}

print_options() {
    echo "Options:"
    echo "    -h               Prints the help message"
    echo "    -e ENVIRONMENT   The environment to start. 'local' is the default option."
    echo "                     The other options are 'prod' and 'test'"
}

while getopts ":e:h" opt; do
    case $opt in
        e)
            env=$OPTARG
            ;;
        h)
            print_description
            echo
            print_options
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_options >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            print_options >&2
            exit 1
            ;;
    esac
done

if [[ -z $env ]]; then
    env="local"
fi

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

# Load get_ip_address function.
source common.sh
export HOST_IP=$(get_ip_address)

# Use a different Nomad port for test environment to avoid interfering
# with the non-test Nomad.
export NOMAD_PORT=4646
if [ $env == "test" ]; then
    export NOMAD_PORT=5646

    # format_nomad_with_env.sh will create distinct test Nomad job
    # specifications to avoid overwriting other job
    # specifications. This is done by appending '.test' to the job
    # specification file names.
    export TEST_POSTFIX="_test"

    # Additional Nomad configuration is necessary to avoid conflicts
    # with non-test Nomad agent.
    export TEST_NOMAD_CONFIG="-config=test_nomad_config.hcl"
else
    # This is only for running Nomad in non-cloud environments so if
    # its not test then we're in local
    env="local"
fi

nomad_dir="$script_directory/nomad_dir$TEST_POSTFIX"
if [ ! -d $nomad_dir ]; then
    mkdir $nomad_dir
fi

# Start Nomad in both server and client mode locally
nomad agent -bind $HOST_IP \
      -data-dir $nomad_dir $TEST_NOMAD_CONFIG \
      -config nomad_client.config \
      -dev \
    &> nomad.logs"$TEST_POSTFIX" &

export NOMAD_ADDR=http://$HOST_IP:$NOMAD_PORT

# While we wait for Nomad to start, make sure the Docker registry has
# an up-to-date Docker image for the workers sub-project. We run a
# local Docker registry because Nomad must pull from a registry, it
# will not use an image which lives on the host machine.

# Make sure the local Docker registry is running.
if [[ -z $(docker ps | grep "image-registry") ]]; then
    docker run -d -p 5000:5000 --restart=always --name image-registry registry:2
fi

# Wait for Nomad to get started.
## A maximum of 1 minute, since it shouldn't take longer than that.
END_SECONDS=$(($SECONDS+60))
nomad_status=$(check_nomad_status)
while [ $nomad_status != "200" ] && [ $SECONDS -lt $END_SECONDS ]; do
    sleep 1
    nomad_status=$(check_nomad_status)
done

echo "Nomad is online. Registering jobs."

./format_nomad_with_env.sh -p workers -e $env
./format_nomad_with_env.sh -p surveyor -e $env

# Register the jobs for dispatching.
for job_spec in $(ls -1 workers/nomad-job-specs/*.nomad$TEST_POSTFIX); do
    echo "Registering $job_spec"
    nomad run $job_spec
done

for job_spec in $(ls -1 foreman/nomad-job-specs/*.nomad$TEST_POSTFIX); do
    echo "Registering $job_spec"
    nomad run $job_spec
done
