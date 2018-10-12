#!/bin/bash -e

cd ~/refinebio

run_on_deploy_box () {
    ssh -o StrictHostKeyChecking=no \
        -i infrastructure/data-refinery-key.pem \
        ubuntu@$DEPLOY_IP_ADDRESS $1
}

# Clear out the old git repository.
run_on_deploy_box "sudo rm -rf refinebio"

# Create file containing local env vars that are needed for deploy.
rm -f env_vars
echo "export CIRCLE_TAG=$CIRCLE_TAG" >> env_vars
echo "export DOCKER_ID=$DOCKER_ID" >> env_vars
echo "export DOCKER_PASSWD=$DOCKER_PASSWD" >> env_vars

# Clear out test data so we don't transfer several GB we don't need
rm -rf workers/test_volume

# # Move the un-git-crypted repo which is on the correct commit to the deploy box.
scp -o StrictHostKeyChecking=no \
    -i infrastructure/data-refinery-key.pem \
    -r . ubuntu@$DEPLOY_IP_ADDRESS:refinebio

run_on_deploy_box "source refinebio/env_vars && sudo echo -e '######\nStarting new deploy for $CIRCLE_TAG\n######' >> /var/log/docker_update.log"
run_on_deploy_box "source refinebio/env_vars && bash refinebio/.circleci/update_docker_img.sh >> /var/log/docker_update.log"
run_on_deploy_box "source refinebio/env_vars && sudo echo -e '######\nStarting new deploy for $CIRCLE_TAG\n######' >> /var/log/deploy.log"
run_on_deploy_box "source refinebio/env_vars && bash refinebio/.circleci/run_terraform.sh >> /var/log/deploy.log"
