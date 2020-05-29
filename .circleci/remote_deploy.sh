#!/bin/bash -e

# This script will deploy the Refinebio system using a dedicated AWS instance.
#   First it will use that instance to build up to date Docker images
#     and push them to Dockerhub.
#   Next it will use terraform to update our infrastructure and restart our services.
#   Finally it cleans up after itself.

# It has been written with the intention of being run from CircleCI as
#   part of our CI/CD process. It therefore assumes that the following
#   environment variables will be set:
#     - DEPLOY_IP_ADDRESS --  The IP address of the instance to run the deploy on.
#     - CIRCLE_TAG -- The tag that was pushed to CircleCI to trigger the deploy.
#         Will be used as the version for the system and the tag for Docker images.
#     - DOCKER_ID -- The username that will be used to log into Dockerhub.
#     - DOCKER_PASSWD -- The password that will be used to log into Dockerhub.
#     - OPENSSL_KEY -- The OpenSSl key which will be used to decrypt the SSH key.
#     - AWS_ACCESS_KEY_ID -- The AWS key id to use when interacting with AWS.
#     - AWS_SECRET_ACCESS_KEY -- The AWS secret key to use when interacting with AWS.


cd ~/refinebio

chmod 600 infrastructure/data-refinery-key.pem

run_on_deploy_box () {
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -i infrastructure/data-refinery-key.pem \
        ubuntu@"${DEPLOY_IP_ADDRESS}" "cd refinebio && $1"
}

# Create file containing local env vars that are needed for deploy.
rm -f env_vars
echo "export CIRCLE_TAG=$CIRCLE_TAG" >> env_vars
echo "export DOCKER_ID=$DOCKER_ID" >> env_vars
echo "export DOCKER_PASSWD=$DOCKER_PASSWD" >> env_vars
echo "export OPENSSL_KEY=\"$OPENSSL_KEY\"" >> env_vars
echo "export AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID" >> env_vars
echo "export AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY" >> env_vars

# And checkout the correct tag.
run_on_deploy_box "git fetch --all"
run_on_deploy_box "git checkout $CIRCLE_TAG"

# Verify that the tag has been signed by a trusted team member.
run_on_deploy_box "bash .circleci/verify_tag.sh $CIRCLE_TAG"

# Copy the necessary environment variables over.
scp -o StrictHostKeyChecking=no \
    -i infrastructure/data-refinery-key.pem \
    -r env_vars ubuntu@"$DEPLOY_IP_ADDRESS":refinebio/env_vars

# Decrypt the secrets in our repo.
run_on_deploy_box "source env_vars && bash .circleci/git_decrypt.sh"

# Output to CircleCI
echo "Building new images"
# Output to the docker update log.
run_on_deploy_box "sudo touch /var/log/docker_update_$CIRCLE_TAG.log"
run_on_deploy_box "sudo chown ubuntu:ubuntu /var/log/docker_update_$CIRCLE_TAG.log"
run_on_deploy_box "source env_vars && echo -e '######\nBuilding new images for $CIRCLE_TAG\n######' 2>&1 | tee -a /var/log/docker_update_$CIRCLE_TAG.log"
run_on_deploy_box "source env_vars && ./.circleci/update_docker_img.sh 2>&1 | tee -a /var/log/docker_update_$CIRCLE_TAG.log"
run_on_deploy_box "source env_vars && echo -e '######\nFinished building new images for $CIRCLE_TAG\n######' 2>&1 | tee -a /var/log/docker_update_$CIRCLE_TAG.log"

# Load docker_img_exists function and $ALL_CCDL_IMAGES
source ~/refinebio/scripts/common.sh

# Circle won't set the branch name for us, so do it ourselves.
branch=$(get_master_or_dev $CIRCLE_TAG)

if [[ "$branch" == "master" ]]; then
    DOCKERHUB_REPO=ccdl
elif [[ "$branch" == "dev" ]]; then
    DOCKERHUB_REPO=ccdlstaging
else
    echo "Why in the world was remote_deploy.sh called from a branch other than dev or master?!?!?"
    exit 1
fi

# It's somehow possible for Docker to sometimes not successfully push
# an image but yet still exit successfully. See:
# https://github.com/AlexsLemonade/refinebio/issues/784
# Since it's not clear how that happened, the safest thing is to add
# an explicit check that the Docker images were successfully updated.
for IMAGE in $ALL_CCDL_IMAGES; do
    image_name="$DOCKERHUB_REPO/dr_$IMAGE"
    if ! docker_img_exists "$image_name" "$CIRCLE_TAG"; then
        echo "Docker image $image_name:$CIRCLE_TAG doesn't exist after running update_docker_img.sh!"
        echo "This is generally caused by a temporary error, please try the 'Rerun job with SSH' button."
        exit 1
    fi
done

# Notify CircleCI that the images have been built.
echo "Finished building new images, running run_terraform.sh."

run_on_deploy_box "sudo touch /var/log/deploy_$CIRCLE_TAG.log"
run_on_deploy_box "sudo chown ubuntu:ubuntu /var/log/deploy_$CIRCLE_TAG.log"
run_on_deploy_box "source env_vars && echo -e '######\nStarting new deploy for $CIRCLE_TAG\n######' >> /var/log/deploy_$CIRCLE_TAG.log 2>&1"
run_on_deploy_box "source env_vars && ./.circleci/run_terraform.sh >> /var/log/deploy_$CIRCLE_TAG.log 2>&1"
run_on_deploy_box "source env_vars && echo -e '######\nDeploying $CIRCLE_TAG finished!\n######' >> /var/log/deploy_$CIRCLE_TAG.log 2>&1"

# Don't leave secrets lying around.
## Clean out any files we've created or moved so git-crypt will relock the repo.
run_on_deploy_box "git clean -f"
run_on_deploy_box "git-crypt lock"
