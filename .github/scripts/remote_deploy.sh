#!/bin/bash -e

# This script will deploy the Refinebio system using a dedicated AWS instance.
#   First it will use that instance to build up to date Docker images
#     and push them to Dockerhub.
#   Next it will use terraform to update our infrastructure and restart our services.
#   Finally it cleans up after itself.

# It has been written with the intention of being run from GitHub Actions as
#   part of our CI/CD process. It therefore assumes that the following
#   environment variables will be set:
#     - DEPLOY_BOX_IP --  The IP address of the instance to run the deploy on.
#     - DEPLOY_TAG -- The tag that was pushed to GitHub to trigger the deploy.
#         Will be used as the version for the system and the tag for Docker images.
#     - DEPLOY_USER -- The user that triggered the deploy (for resource tagging).
#     - BRANCH -- The deploy branch (master|dev), resolved upstream by the
#         determine_branch job.
#     - DOCKER_USERNAME -- The username that will be used to log into Dockerhub.
#     - DOCKER_PASSWORD -- The password that will be used to log into Dockerhub.
#     - AWS_ACCESS_KEY_ID -- The AWS key id to use when interacting with AWS.
#     - AWS_SECRET_ACCESS_KEY -- The AWS secret key to use when interacting with AWS.

echo "$DEPLOY_BOX_SSH_PRIVATE_KEY" >infrastructure/data-refinery-key.pem
chmod 600 infrastructure/data-refinery-key.pem

run_on_deploy_box() {
    # shellcheck disable=SC2029
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -i infrastructure/data-refinery-key.pem \
        ubuntu@"${DEPLOY_BOX_IP}" \
        "cd refinebio && $1"
}

# Resolve env-specific values from the branch determined upstream. These are
# written into env_vars (below) so the deploy box's rbio commands can read them.
case "$BRANCH" in
    master)
        DEPLOY_ENV=prod
        DOCKERHUB_REPO=ccdl
        BATCH_USE_ON_DEMAND_INSTANCES=false
        ;;
    dev)
        DEPLOY_ENV=staging
        DOCKERHUB_REPO=ccdlstaging
        BATCH_USE_ON_DEMAND_INSTANCES=true
        ;;
    *)
        echo "remote_deploy.sh: unexpected BRANCH=$BRANCH (expected master or dev)" >&2
        exit 1
        ;;
esac

# Create file containing local env vars that are needed for deploy.
# SYSTEM_VERSION mirrors DEPLOY_TAG — docker-bake.hcl reads it as the image tag.
# PUSH_CACHE=true tells the bake config to also publish cache layers to the
# registry so subsequent deploy builds are warm.
rm -f env_vars
cat >>env_vars <<EOF
export DEPLOY_TAG='$DEPLOY_TAG'
export DEPLOY_USER='$DEPLOY_USER'
export DEPLOY_ENV='$DEPLOY_ENV'
export SYSTEM_VERSION='$DEPLOY_TAG'
export DOCKERHUB_REPO='$DOCKERHUB_REPO'
export PUSH_CACHE=true
export BATCH_USE_ON_DEMAND_INSTANCES='$BATCH_USE_ON_DEMAND_INSTANCES'
export DOCKER_USERNAME='$DOCKER_USERNAME'
export DOCKER_PASSWORD='$DOCKER_PASSWORD'
export AWS_ACCESS_KEY_ID='$AWS_ACCESS_KEY_ID'
export AWS_SECRET_ACCESS_KEY='$AWS_SECRET_ACCESS_KEY'
export TF_VAR_database_password='$DATABASE_PASSWORD'
export TF_VAR_django_secret_key='$DJANGO_SECRET_KEY'
export TF_VAR_sentry_dsn='$SENTRY_DSN'
export TF_VAR_slack_webhook_url='$SLACK_WEBHOOK_URL'
export TF_VAR_ssh_public_key='$SSH_PUBLIC_KEY'
EOF

# And checkout the correct tag.
run_on_deploy_box "git fetch --all"
run_on_deploy_box "git checkout $DEPLOY_TAG"

# Verify that the tag has been signed by a trusted team member.
run_on_deploy_box "bash .github/scripts/verify_tag.sh $DEPLOY_TAG"

# Copy the necessary environment variables over.
scp -o StrictHostKeyChecking=no \
    -i infrastructure/data-refinery-key.pem \
    -r env_vars ubuntu@"$DEPLOY_BOX_IP":refinebio/env_vars

# Along with the ssh key iself, which the deploy script will use.
scp -o StrictHostKeyChecking=no \
    -i infrastructure/data-refinery-key.pem \
    -r infrastructure/data-refinery-key.pem ubuntu@"$DEPLOY_BOX_IP":refinebio/infrastructure/data-refinery-key.pem

echo "Building new images"
run_on_deploy_box "sudo touch /var/log/docker_update_$DEPLOY_TAG.log"
run_on_deploy_box "sudo chown ubuntu:ubuntu /var/log/docker_update_$DEPLOY_TAG.log"
run_on_deploy_box ". env_vars && echo -e '######\nBuilding new images for $DEPLOY_TAG\n######' 2>&1 | tee -a /var/log/docker_update_$DEPLOY_TAG.log"
# DOCKER_IO_USERNAME / DOCKER_IO_PASSWORD are expected to be set in the deploy box's environment. They are probably legacy.
run_on_deploy_box "docker login -u \"\$DOCKER_IO_USERNAME\" -p \"\$DOCKER_IO_PASSWORD\""
run_on_deploy_box ". env_vars && ./bin/rbio build --push deploy 2>&1 | tee -a /var/log/docker_update_$DEPLOY_TAG.log"
run_on_deploy_box ". env_vars && echo -e '######\nFinished building new images for $DEPLOY_TAG\n######' 2>&1 | tee -a /var/log/docker_update_$DEPLOY_TAG.log"

# Docker sometimes exits 0 from `bake --push` even when the push failed (see
# https://github.com/AlexsLemonade/refinebio/issues/784). Verify each image
# is actually on Dockerhub before proceeding to deploy.
docker_image_exists() {
    TOKEN=$(curl -s -H "Content-Type: application/json" -X POST \
        -d '{"username": "'"${DOCKER_USERNAME}"'", "password": "'"${DOCKER_PASSWORD}"'"}' \
        https://hub.docker.com/v2/users/login/ | jq -r .token)
    EXISTS=$(curl -s -H "Authorization: JWT ${TOKEN}" \
        "https://hub.docker.com/v2/repositories/$1/tags/?page_size=10000" |
        jq -r "[.results | .[] | .name == \"$2\"] | any" 2>/dev/null)
    test -n "$EXISTS" -a "$EXISTS" = true
}

ALL_IMAGES="base api_base api foreman smasher compendia illumina affymetrix salmon transcriptome no_op downloaders"

for image in $ALL_IMAGES; do
    image_name="$DOCKERHUB_REPO/dr_$image"
    if ! docker_image_exists "$image_name" "$DEPLOY_TAG"; then
        echo "Docker image $image_name:$DEPLOY_TAG doesn't exist after building!"
        echo "This is generally caused by a temporary error, please try the 'Rerun workflow' button."
        exit 1
    fi
done

echo "Finished building new images, running rbio deploy:up."

run_on_deploy_box "sudo touch /var/log/deploy_$DEPLOY_TAG.log"
run_on_deploy_box "sudo chown ubuntu:ubuntu /var/log/deploy_$DEPLOY_TAG.log"
run_on_deploy_box ". env_vars && echo -e '######\nStarting new deploy for $DEPLOY_TAG\n######' >> /var/log/deploy_$DEPLOY_TAG.log 2>&1"
run_on_deploy_box "sudo .github/scripts/update_ca_certificates.sh >> /var/log/deploy_$DEPLOY_TAG.log 2>&1"

# This should never be logged to Github Actions in case it exposes any secrets as terraform output.
run_on_deploy_box ". env_vars && ./bin/rbio deploy:up >> /var/log/deploy_$DEPLOY_TAG.log 2>&1"

run_on_deploy_box ". env_vars && echo -e '######\nDeploying $DEPLOY_TAG finished!\n######' >> /var/log/deploy_$DEPLOY_TAG.log 2>&1"

.github/scripts/slackpost_deploy.sh robots deploybot

# Temporarily disabled due to surveyor pipeline, only API + smasher are currently operational.
# if [[ "$BRANCH" == "dev" ]]; then
#     run_on_deploy_box ". env_vars && ./foreman/run_end_to_end_tests.sh"
#     .github/scripts/slackpost_end_to_end.sh robots deploybot
# fi

# Don't leave secrets lying around.
run_on_deploy_box "git clean -f"
