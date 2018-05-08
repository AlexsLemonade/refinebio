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
            echo "Starts Nomad and registers the jobs with it. This involves re-building the "
            echo "Docker images and running format_nomad_with_env.sh to format the Nomad job "
            echo "specifications for the correct environment."
            echo '- "dev" is the default enviroment, use -e to specify "prod" or "test".'
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

# Load get_ip_address and docker_img_exists functions.
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
    env="dev"
fi

nomad_dir="$script_directory/nomad_dir$TEST_POSTFIX"
if [ ! -d $nomad_dir ]; then
    mkdir $nomad_dir
fi

echo "Pulling/rebuilding Docker images while waiting for Nomad to come online."

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

branch_name=$(git rev-parse --symbolic-full-name --abbrev-ref HEAD)

# We want to check if a test image has been built for this branch. If
# it has we should use that rather than building it slowly. In the
# case of affymetrix though, the image takes too long to build so if
# they haven't built it, just use the latest prod image.

worker_images=(affymetrix salmon transcriptome no_op downloaders)

# Build, tag, and push an image for the workers to the local registry.
for image in ${worker_images[*]}; do
    image_name=ccdl/dr_$image
    if docker_img_exists $image_name $branch_name ; then
        docker pull $image_name:$branch_name
        # Don't push the affymetrix image because it's huge and will
        # make the tests take forever.
        if [ $image != "ccdl/dr_affymetrix" ]; then
            docker push localhost:5000/"$image_name"
        fi
    elif [ $image_name == "ccdl/dr_affymetrix" ]; then
        docker pull $image_name:latest
    else
        echo ""
        echo "Rebuilding the $image_name image."
        docker build -t "$image_name" -f workers/dockerfiles/Dockerfile.$image .
        docker tag "$image_name" localhost:5000/"$image_name"
        docker push localhost:5000/"$image_name"
    fi
done

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
for job_spec in $(ls -1 workers/nomad-job-specs | grep "\.nomad$TEST_POSTFIX$"); do
    nomad run workers/nomad-job-specs/"$job_spec"
done

# There's only one foreman image, so no need to loop.
nomad run foreman/nomad-job-specs/surveyor.nomad"$TEST_POSTFIX"
