#!/bin/bash -e

# This script must be run from /infrastructure!

# Make our IP address known to terraform.
source ../common.sh
export TF_VAR_host_ip=`wget -qO- http://ipecho.net/plain ; echo`

# Install Nomad if not already installed.
if [[ $(nomad -version |& grep command) != "" ]]; then
    ../install_nomad.sh
fi

# Copy ingress config to top level so it can be applied.
cp deploy/ci_ingress.tf .

# Open up ingress to AWS for Circle, stop jobs, migrate DB.
echo "Applying ingress.."
terraform apply -auto-approve

# Find address of Nomad server.
export NOMAD_LEAD_SERVER_IP=`terraform output nomad_server_1_ip`
export NOMAD_ADDR=http://$NOMAD_LEAD_SERVER_IP:4646

# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  $NOMAD_ADDR/v1/status/leader)
}


# Wait for Nomad to get started in case the server just went up for
# the first time.
echo "Confirming Nomad cluster.."
start_time=$(date +%s)
diff=0
nomad_status=$(check_nomad_status)
while [[ $diff < 300 && $nomad_status != "200" ]]; do
    sleep 1
    nomad_status=$(check_nomad_status)
    let "diff = $(date +%s) - $start_time"
done

# Kill Base Nomad Jobs so no new jobs can be queued.
echo "Killing base jobs.."
if [[ $(nomad status) != "No running jobs" ]]; then
    for job in $(nomad status | grep running | awk {'print $1'} || grep --invert-match /)
    do
        nomad stop $job
    done
fi

# Kill parameterized Nomad Jobs so no jobs will be running when we
# apply migrations.
echo "Killing dispatch jobs.."
if [[ $(nomad status) != "No running jobs" ]]; then
    for job in $(nomad status | awk {'print $1'} || grep /)
    do
        # Skip the header row for jobs.
        if [ $job != "ID" ]; then
            nomad stop $job
        fi
    done
fi

# Declare our production environment variables. These should be in
# parity with the env files such as workers/environments/dev.
# Note:
# I have two options here. Make this the source of truth for these
# variables or have terraform spit every single one of them out as
# outputs.
export REGION=`terraform output region`
export USER=`terraform output user`
export STAGE=`terraform output stage`
export AWS_ACCESS_KEY_ID=`terraform output aws_access_key_id`
export AWS_SECRET_ACCESS_KEY=`terraform output aws_secret_access_key`
export DJANGO_DEBUG=False
export DATABASE_NAME=`terraform output database_name`
export DATABASE_HOST=`terraform output database_host`
export DATABASE_USER=drpostgresuser
# export DATABASE_PASSWORD=drpostgrespassword
export DATABASE_PASSWORD=629f6R4eUgNmBzJf6Qthaw5wzaKqT9UA
export DATABASE_PORT=5432
export DATABASE_TIMEOUT=30
export DJANGO_SECRET_KEY=NtG1bxZU115GThwrLuAJe0PhTVN9hJ4P
export RUNNING_IN_CLOUD=True
export USE_S3=True
export S3_BUCKET_NAME=data-refinery # Should be TF output?
export LOCAL_ROOT_DIR=/home/user/data_store
export RAW_PREFIX=raw
export TEMP_PREFIX=temp
export PROCESSED_PREFIX=processed
export WORKERS_DOCKER_IMAGE=miserlou/dr_worker:2
export FOREMAN_DOCKER_IMAGE=miserlou/dr_foreman:3

# Big hack here. Figure something better out tomorrow.
echo "REGION=`terraform output region`" > prod_env
echo "USER=`terraform output user`" >> prod_env
echo "STAGE=`terraform output stage`" >> prod_env
echo "AWS_ACCESS_KEY_ID=`terraform output aws_access_key_id`" >> prod_env
echo "AWS_SECRET_ACCESS_KEY=`terraform output aws_secret_access_key`" >> prod_env
echo "DJANGO_DEBUG=False" >> prod_env
echo "DATABASE_NAME=`terraform output database_name`" >> prod_env
echo "DATABASE_HOST=`terraform output database_host`" >> prod_env
echo "DATABASE_USER=drpostgresuser" >> prod_env
# Not sure what's going on with this here password, but it doesn't seem to work.
echo "DATABASE_PASSWORD=629f6R4eUgNmBzJf6Qthaw5wzaKqT9UA" >> prod_env
echo "DATABASE_PORT=5432" >> prod_env
echo "DATABASE_TIMEOUT=30" >> prod_env
echo "DJANGO_SECRET_KEY=NtG1bxZU115GThwrLuAJe0PhTVN9hJ4P" >> prod_env
echo "RUNNING_IN_CLOUD=True" >> prod_env
echo "USE_S3=True" >> prod_env
echo "S3_BUCKET_NAME=data-refinery # Should be TF output?" >> prod_env
echo "LOCAL_ROOT_DIR=/home/user/data_store" >> prod_env
echo "RAW_PREFIX=raw" >> prod_env
echo "TEMP_PREFIX=temp" >> prod_env
echo "PROCESSED_PREFIX=processed" >> prod_env
echo "WORKERS_DOCKER_IMAGE=miserlou/dr_worker:2" >> prod_env
echo "FOREMAN_DOCKER_IMAGE=miserlou/dr_foreman:3" >> prod_env

# Create directory for migration files.
echo "Migrating.."
mkdir -p migrations;

# Get an image to run the migrations with.
docker pull miserlou/dr_foreman:3;

# Make the migration files.
docker run \
       --volume migrations \
       --env-file prod_env \
       miserlou/dr_foreman:3 makemigrations;

# Migrate auth.
docker run \
       --volume migrations \
       --env-file prod_env \
       miserlou/dr_foreman:3 migrate auth;

# Apply general migrations
docker run \
       --volume migrations \
       --env-file prod_env \
       miserlou/dr_foreman:3 migrate;

../workers/format_nomad_with_env.sh -e prod -o $(pwd)/nomad-job-specs/
../foreman/format_nomad_with_env.sh -e prod -o $(pwd)/nomad-job-specs/

# Re-register Nomad jobs.
echo "Registering new job specifications.."
nomad_job_specs=nomad-job-specs/*
for nomad_job_spec in $nomad_job_specs; do
    echo "registering $nomad_job_spec"
    nomad run -address http://$NOMAD_LEAD_SERVER_IP:4646 $nomad_job_spec
done

# Remove the ingress config so the next `terraform apply` will remove
# access for Circle.
# rm ci_ingress.tf
