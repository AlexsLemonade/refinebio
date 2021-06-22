#!/bin/sh

# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

./scripts/prepare_image.sh -i smasher
image_name="ccdlstaging/dr_smasher"

volume_directory="$script_directory/volume"

. ./scripts/common.sh
DB_HOST_IP=$(get_docker_db_ip_address)

AWS_ACCESS_KEY_ID="$(~/bin/aws configure get default.aws_access_key_id)"
export AWS_ACCESS_KEY_ID
AWS_SECRET_ACCESS_KEY="$(~/bin/aws configure get default.aws_secret_access_key)"
export AWS_SECRET_ACCESS_KEY

docker run \
       -it \
       -m 500m \
       --add-host=database:"$DB_HOST_IP" \
       --env-file workers/environments/local \
       --env AWS_ACCESS_KEY_ID \
       --env AWS_SECRET_ACCESS_KEY \
       --entrypoint ./manage.py \
       --volume "$volume_directory":/home/user/data_store \
       --link drdb:postgres \
       "$image_name" run_janitor
