#!/bin/bash

# This script runs the Data Refinery Workers locally.
# What that actually means is that it runs Nomad and registers the
# Nomad jobs for the Data Refinery Workers. Nomad is listening for
# messages telling it to run Processor and Downloader jobs. When it
# receives them, it will run the jobs within new Docker containers.

while getopts ":p:e:o:h" opt; do
    case $opt in
        e)
            env=$OPTARG
            ;;
        h)
            echo "Formats Nomad Job Specifications with the specified environment overlaid "
            echo "onto the current environment."
            echo '-p specifies the project to format. Valid values are "api", workers" or "foreman".'
            echo '- "dev" is the default enviroment, use -e to specify "prod" or "test".'
            echo '- the project directory will be used as the default output directory, use -o to specify'
            echo '      an absolute path to a directory (trailing / must be included).'
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [[ -z $env ]]; then
    # XXX: for now dev==local and prod==cloud. This works because we
    # don't have a true prod environment yet so using prod for cloud
    # development is okay, but we definitely need to address
    # https://github.com/AlexsLemonade/refinebio/issues/199 before we
    # create an actual prod environment.
    env="dev"
fi

# Figure out the right location to put the nomad directory.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/workers/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

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
    # its not test then we're in dev (XXX: change this to local.)
    env = "dev"
fi


nomad_dir="$script_directory/nomad_dir$TEST_POSTFIX"
if [ ! -d $nomad_dir ]; then
    mkdir $nomad_dir
fi

echo "Rebuilding Docker image while waiting for Nomad to come online."

# Start Nomad in both server and client mode locally
nomad agent -bind $HOST_IP \
      -data-dir $nomad_dir $TEST_NOMAD_CONFIG \
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

# Build, tag, and push an image for the workers to the local registry.
docker build -t dr_worker"$TEST_POSTFIX" -f workers/Dockerfile .
docker tag dr_worker"$TEST_POSTFIX" localhost:5000/dr_worker"$TEST_POSTFIX"
docker push localhost:5000/dr_worker"$TEST_POSTFIX"

# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  $NOMAD_ADDR/v1/status/leader)
}

# Wait for Nomad to get started.
nomad_status=$(check_nomad_status)
while [ $nomad_status != "200" ]; do
    sleep 1
    nomad_status=$(check_nomad_status)
done

echo "Nomad is online. Registering jobs."

./format_nomad_with_env.sh -p workers -e $env
./format_nomad_with_env.sh -p foreman -e $env

# Register the jobs for dispatching.
nomad run workers/downloader.nomad"$TEST_POSTFIX"
nomad run workers/processor.nomad"$TEST_POSTFIX"
nomad run foreman/surveyor.nomad"$TEST_POSTFIX"
