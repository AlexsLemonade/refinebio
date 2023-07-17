#!/bin/sh

# The directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# However, in order to give Docker access to all the code we have to
# move up a level.
cd ..

./scripts/prepare_image.sh -i smasher

. ./scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)

AWS_ACCESS_KEY_ID="$(~/bin/aws configure get default.aws_access_key_id)"
export AWS_ACCESS_KEY_ID
AWS_SECRET_ACCESS_KEY="$(~/bin/aws configure get default.aws_secret_access_key)"
export AWS_SECRET_ACCESS_KEY

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --entrypoint ./manage.py \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    --env-file workers/environments/local \
    --interactive \
    --link drdb:postgres \
    --memory 500m \
    --tty \
    --volume "$script_directory/volume":/home/user/data_store \
    "$DOCKERHUB_REPO/dr_smasher" \
    run_janitor
