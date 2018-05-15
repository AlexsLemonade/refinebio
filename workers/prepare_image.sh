#!/bin/bash

while getopts ":i:l:h" opt; do
    case $opt in
        i)
            image=$OPTARG
            ;;
        p)
            pull="True"
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

source common.sh

# We want to check if a test image has been built for this branch. If
# it has we should use that rather than building it slowly.
image_name=ccdl/dr_$image
if [[ "$(docker_img_exists $image_name $branch_name)" ]] ; then
    docker pull $image_name:$branch_name
elif [[ ! -z $pull ]]; then
    docker pull $image_name
else
    echo ""
    echo "Rebuilding the $image_name image."
    docker build -t "$image_name" -f workers/dockerfiles/Dockerfile.$image .
fi
