#!/bin/bash

while getopts "phi:" opt; do
    case $opt in
        i)
            image=$OPTARG
            ;;
        p)
            pull="True"
            ;;
        h)
            echo "Prepares an workers image specified by -i."
            echo "First, it checks to see if there is an existing image on Dockerhub"
            echo "which is tagged with the current branch name. If there is it pulls that."
            echo "Otherwise it will build the image locally, unless instructed"
            echo "to pull the latest version of the image from Dockerhub with -p."
            exit 0
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
image_name="ccdl/dr_$image"
if [[ "$(docker_img_exists $image_name $branch_name)" ]] ; then
    docker pull $image_name:$branch_name
elif [[ ! -z $pull ]]; then
    docker pull $image_name
else
    echo ""
    echo "Rebuilding the $image_name image."
    docker build -t "$image_name" -f workers/dockerfiles/Dockerfile.$image .
fi
