#!/bin/bash -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

print_description() {
    # shellcheck disable=SC2016
    echo 'This script can be used to deploy and update a `refine.bio` instance stack.'
    echo 'It will create all of the AWS infrasctructure (roles/instances/db/network/etc),'
    echo 'open an ingress, kill all running Nomad jobs, perform a database migration,'
    echo 're-define and re-register Nomad job specifications, and finally close the'
    echo 'ingress. This can be run from a CI/CD machine or a local dev box.'
}

print_options() {
    echo 'This script accepts the following arguments: -e, -v, -u, -r, and -h.'
    echo '-h prints this help message and exits.'
    echo '-e specifies the environment you would like to deploy to and is not optional. Its valid values are:'
    echo '   "-e prod" will deploy the production stack. This should only be used from a CD machine.'
    echo '   "-e staging" will deploy the staging stack. This should only be used from a CD machine.'
    echo '   "-e dev" will deploy a dev stack which is appropriate for a single developer to use to test.'
    echo '-d May be used to override the Dockerhub repo where the images will be pulled from.'
    echo '   This may also be specified by setting the TF_VAR_dockerhub_repo environment variable.'
    echo '   If unset, defaults to "ccdlstaging" if the version contains "-dev" and "ccdl" otherwise.'
    echo '   for dev and staging environments and "ccdl" for prod.'
    echo '   This option is useful for testing code changes. Images with the code to be tested can be pushed'
    echo '   to your private Dockerhub repo and then the system will find them.'
    echo '-v specifies the version of the system which is being deployed and is not optional.'
    echo "-u specifies the username of the deployer. Should be the developer's name in development stacks."
    echo '   This option may be omitted, in which case the TF_VAR_user variable MUST be set instead.'
    echo '-r specifies the AWS region to deploy the stack to. Defaults to us-east-1.'
}

while getopts ":e:d:v:u:r:h" opt; do
    case $opt in
    e)
        export env=$OPTARG
        export TF_VAR_stage=$OPTARG
        ;;
    d)
        export TF_VAR_dockerhub_repo=$OPTARG
        ;;
    v)
        export SYSTEM_VERSION=$OPTARG
        ;;
    u)
        export TF_VAR_user=$OPTARG
        ;;
    r)
        export TF_VAR_region=$OPTARG
        ;;
    h)
        print_description
        echo
        print_options
        exit 0
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        print_options >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        print_options >&2
        exit 1
        ;;
    esac
done

if [[ $env != "dev" && $env != "staging" && $env != "prod" ]]; then
    echo 'Error: must specify environment as either "dev", "staging", or "prod" with -e.'
    exit 1
fi

if [[ -z $TF_VAR_user ]]; then
    echo 'Error: must specify the username by either providing the -u argument or setting TF_VAR_user.'
    exit 1
fi

if [[ -z $SYSTEM_VERSION ]]; then
    echo 'Error: must specify the system version with -v.'
    exit 1
fi

if [[ -z $TF_VAR_dockerhub_repo ]]; then
    if [[ $SYSTEM_VERSION == *"-dev" ]]; then
        export TF_VAR_dockerhub_repo=ccdlstaging
    else
        export TF_VAR_dockerhub_repo=ccdl
    fi
fi

if [[ -z $TF_VAR_region ]]; then
    TF_VAR_region=us-east-1
fi

# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    curl --write-out "%{http_code}" \
      --silent \
      --output /dev/null \
      "$NOMAD_ADDR/v1/status/leader"
}

# We have terraform output environment variables via a single output
# variable, which we then read in as json using the command line tool
# `jq`, so that we can use them via bash.
format_environment_variables () {
  json_env_vars=$(terraform output -json environment_variables | jq -c '.value[]')
  for row in $json_env_vars; do
      env_var_assignment=$(echo "$row" | jq -r ".name")=$(echo "$row" | jq -r ".value")
      export "${env_var_assignment?}"
      echo "$env_var_assignment" >> prod_env
  done
}

# Load $ALL_CCDL_IMAGES and helper functions
source ../scripts/common.sh
# Make our IP address known to terraform.
TF_VAR_host_ip="$(dig +short myip.opendns.com @resolver1.opendns.com)"
export TF_VAR_host_ip

for IMAGE in $ALL_CCDL_IMAGES; do
    # For each image we need to set the env var that is used by our
    # scripts and the env var that gets picked up by terraform because
    # it is preceeded with TF_VAR.
    IMAGE_UPPER="$(echo "$IMAGE" | tr '[:lower:]' '[:upper:]')"
    export "${IMAGE_UPPER}_DOCKER_IMAGE=dr_$IMAGE:$SYSTEM_VERSION"
    export "TF_VAR_${IMAGE}_docker_image=dr_$IMAGE:$SYSTEM_VERSION"
done

# Copy ingress config to top level so it can be applied.
cp deploy/ci_ingress.tf .

# Check if a new ccdl-ubuntu ami will be needed for this region
if [[ $(aws ec2 describe-images \
            --region $TF_VAR_region --owners 589864003899 \
            --filters 'Name=name,Values=ccdl-ubuntu-18.04-*' \
            --query 'length(Images)') \
            -eq 0 ]]; then
    echo "No ccdl-ubuntu-18.04 AMI found for this region, creating a new one"

    # Find most recent ccdl-ubuntu ami from us-east-1
    template_ami_id=$(aws ec2 describe-images \
                          --region us-east-1 --owners 589864003899 \
                          --filters 'Name=name,Values=ccdl-ubuntu-18.04-*' \
                          --query 'sort_by(Images,&CreationDate)[-1].ImageId' \
                          --output text)
    
    # Make a copy into this region
    new_ami_name="ccdl-ubuntu-18.04-$(date "+%Y-%m-%dT%H.%M.%S")"
    new_ami_id=$(aws ec2 copy-image --source-image-id $template_ami_id --source-region us-east-1 --region $TF_VAR_region --name $new_ami_name --output text)
    echo "Created new AMI for $TF_VAR_region"
    echo "    name: $new_ami_name"
    echo "    id:   $new_ami_id"
fi

# Always init terraform first, especially since we're using a remote backend.
./init_terraform.sh

if [[ ! -f terraform.tfstate ]]; then
    ran_init_build=true
    echo "No terraform state file found, applying initial terraform deployment."

    # These files are inputs but are created by format_nomad_with_env.sh
    # based on outputs from terraform. Kinda a Catch 22, but we can
    # get around it by providing dummy files to get bootstrapped.
    touch api-configuration/environment
    touch foreman-configuration/environment

    # Output the plan for debugging deployments later.
    # Until terraform plan supports -var-file the plan is wrong.
    # terraform plan

    if [[ -n $CIRCLE_BUILD_NUM ]]; then
        # Make sure we can't expose secrets in circleci
        terraform apply -var-file="environments/$env.tfvars" -auto-approve > /dev/null
    else
        terraform apply -var-file="environments/$env.tfvars" -auto-approve
    fi
fi

# We have to do this once before the initial deploy..
rm -f prod_env
format_environment_variables

../scripts/format_nomad_with_env.sh -p api -e "$env" -o "$(pwd)/api-configuration/"
../scripts/format_nomad_with_env.sh -p foreman -e "$env" -o "$(pwd)/foreman-configuration/"

if [[ -z $ran_init_build ]]; then
    # Open up ingress to AWS for Circle, stop jobs, migrate DB.
    echo "Deploying with ingress.."

    # Output the plan for debugging deployments later.
    # Until terraform plan supports -var-file the plan is wrong.
    # terraform plan

    if [[ -n $CIRCLE_BUILD_NUM ]]; then
        # Make sure we can't expose secrets in circleci
        terraform apply -var-file="environments/$env.tfvars" -auto-approve > /dev/null
    else
        terraform apply -var-file="environments/$env.tfvars" -auto-approve
    fi
fi

# Find address of Nomad server.
NOMAD_LEAD_SERVER_IP="$(terraform output nomad_server_1_ip)"
export NOMAD_LEAD_SERVER_IP

export NOMAD_ADDR=http://$NOMAD_LEAD_SERVER_IP:4646

# Wait for Nomad to get started in case the server just went up for
# the first time.

set +e # curl fails if the nomad server isn't up

echo "Confirming Nomad cluster.."
start_time=$(date +%s)
diff=0
nomad_status=$(check_nomad_status)
while [ "$diff" -lt "900" ] && [ "$nomad_status" != "200" ]; do
    sleep 1
    nomad_status=$(check_nomad_status)
    (( diff = $(date +%s) - start_time ))
done

if [[ $nomad_status != "200" ]]; then
    echo "Nomad didn't start, aborting deploy."
    echo "Either the timeout needs to be raised or there's something else broken."
    exit 1
fi

set -e # Turn errors back on after we confirm that the nomad server is up

# Kill Base Nomad Jobs so no new jobs can be queued.
echo "Killing base jobs.. (this takes a while..)"
if [[ "$(nomad status)" != "No running jobs" ]]; then
    for job in $(nomad status | grep running | awk '{print $1}' | grep --invert-match /)
    do
        # '|| true' so that if a job is garbage collected before we can remove it the error
        # doesn't interrupt our deploy.
        nomad stop -purge -detach "$job" > /dev/null || true &
    done
fi

# Wait to make sure that all base jobs are killed so no new jobs can
# be queued while we kill the parameterized Nomad jobs.
# shellcheck disable=SC2046
wait $(jobs -p)

# Kill parameterized Nomad Jobs so no jobs will be running when we
# apply migrations.
echo "Killing dispatch jobs.. (this also takes a while..)"
if [[ "$(nomad status)" != "No running jobs" ]]; then
    counter=0
    for job in $(nomad status | awk '{print $1}' | grep /)
    do
        # Skip the header row for jobs.
        if [ "$job" != "ID" ]; then
            # '|| true' so that if a job is garbage collected before we can remove it the error
            # doesn't interrupt our deploy.
            nomad stop -purge -detach "$job" > /dev/null || true &
            counter=$((counter+1))
        fi

        # Wait for all the jobs to stop every 100 so we don't knock
        # over the deploy box if there are 1000's.
        if [[ "$counter" -gt 100 ]]; then
            # shellcheck disable=SC2046
            wait $(jobs -p)
            counter=0
        fi
    done
fi

# Wait for any remaining jobs to all die.
# shellcheck disable=SC2046
wait $(jobs -p)

# Make sure that prod_env is empty since we are only appending to it.
# prod_env is a temporary file we use to pass environment variables to
# `docker run` commands when running migrations.
rm -f prod_env

# (cont'd) ..and once again after the update when this is re-run.
format_environment_variables

# Get an image to run the migrations with.
docker pull "$DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE"

# Test that the pg_bouncer instance is up. 15 minutes should be more than enough.
start_time=$(date +%s)
diff=0
until pg_isready -d "$DATABASE_NAME" -h "$DATABASE_PUBLIC_HOST" -p "$DATABASE_PORT" -U "$DATABASE_USER" &> /dev/null || [ "$diff" -gt "900" ]
do
    echo "Waiting for the pg_bouncer instance to come online ..."
    sleep 10
    (( diff = $(date +%s) - start_time ))
done

if ! pg_isready -d "$DATABASE_NAME" -h "$DATABASE_PUBLIC_HOST" -p "$DATABASE_PORT" -U "$DATABASE_USER" &> /dev/null; then
    echo "pg_bouncer instance failed to come up after 15 minutes."
    exit 1
fi

# Migrate auth.
docker run \
       --env-file prod_env \
       --env RUNNING_IN_CLOUD=False \
       --env DATABASE_HOST="$DATABASE_PUBLIC_HOST" \
       "$DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE" python3 manage.py migrate auth

# Apply general migrations.
docker run \
       --env-file prod_env \
       --env RUNNING_IN_CLOUD=False \
       --env DATABASE_HOST="$DATABASE_PUBLIC_HOST" \
       "$DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE" python3 manage.py migrate

# Create the cache table if it does not already exist.
docker run \
       --env-file prod_env \
       --env RUNNING_IN_CLOUD=False \
       --env DATABASE_HOST="$DATABASE_PUBLIC_HOST" \
       "$DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE" python3 manage.py createcachetable

# Make sure to clear out any old nomad job specifications since we
# will register everything in this directory.
if [ -e nomad-job-specs ]; then
  rm -r nomad-job-specs
fi

# Template the environment variables for production into the Nomad Job
# specs and API confs.
mkdir -p nomad-job-specs
../scripts/format_nomad_with_env.sh -p workers -e "$env" -o "$(pwd)/nomad-job-specs"
../scripts/format_nomad_with_env.sh -p surveyor -e "$env" -o "$(pwd)/nomad-job-specs"

# API and foreman aren't run as nomad jobs, but the templater still works.
../scripts/format_nomad_with_env.sh -p foreman -e "$env" -o "$(pwd)/foreman-configuration"
../scripts/format_nomad_with_env.sh -p api -e "$env" -o "$(pwd)/api-configuration/"

# Re-register Nomad jobs (skip those that end in .tpl)
echo "Registering new job specifications.."
for nomad_job_spec in nomad-job-specs/*.nomad; do
    nomad run "$nomad_job_spec" &
done
echo "Job registrations have been fired off."

# Prepare the client instance user data script for the nomad client instances.
# The `prepare-client-instance-user-data.sh` script overwrites
# `client-instance-user-data.tpl.sh`, so we have to back it up first.
if [ ! -f client-instance-user-data.tpl.sh.bak ]; then
    mv nomad-configuration/client-instance-user-data.tpl.sh nomad-configuration/client-instance-user-data.tpl.sh.bak
fi

./nomad-configuration/prepare-client-instance-user-data.sh
terraform taint aws_spot_fleet_request.cheap_ram

# Ensure the latest image version is being used for the Foreman
terraform taint aws_instance.foreman_server_1

# Remove the ingress config so the next `terraform apply` will remove
# access for Circle.
echo "Removing ingress.."
rm ci_ingress.tf

if [[ -n $CIRCLE_BUILD_NUM ]]; then
    # Make sure we can't expose secrets in circleci
    terraform apply -var-file="environments/$env.tfvars" -auto-approve > /dev/null
else
    terraform apply -var-file="environments/$env.tfvars" -auto-approve
fi

# Don't leave secrets lying around!
rm -f prod_env
# The tarball at the end of client-instance-user-data.tpl.sh has secrets, so restore from backup
mv nomad-configuration/client-instance-user-data.tpl.sh.bak nomad-configuration/client-instance-user-data.tpl.sh

# We try to avoid rebuilding the API server because we can only run certbot
# 5 times a week. Therefore we pull the newest image and restart the API
# this way rather than by tainting the server like we do for foreman.
chmod 600 data-refinery-key.pem
API_IP_ADDRESS=$(terraform output -json api_server_1_ip | jq -c '.value' | tr -d '"')

# To check to see if the docker container needs to be stopped before
# it can be started, grep for the name of the container. However if
# it's not found then grep will return a non-zero exit code so in that
# case return an empty string.
container_running=$(ssh -o StrictHostKeyChecking=no \
                        -o ServerAliveInterval=15 \
                        -o ConnectTimeout=5 \
                        -i data-refinery-key.pem \
                        "ubuntu@$API_IP_ADDRESS"  "docker ps -a" | grep dr_api || echo "")

# If $container_running is empty, then it's because the container isn't running.
# If the container isn't running, then it's because the instance is spinning up.
# The container will be started by the API's init script, so no need to do anything more.

# However if $container_running isn't empty then we need to stop and restart it.
if [[ -n $container_running ]]; then
    echo "Restarting API with latest image."

    # shellcheck disable=SC2029
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        "ubuntu@$API_IP_ADDRESS" "docker pull $DOCKERHUB_REPO/$API_DOCKER_IMAGE"

    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        "ubuntu@$API_IP_ADDRESS" "docker rm -f dr_api"

    scp -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        api-configuration/environment "ubuntu@$API_IP_ADDRESS:/home/ubuntu/environment"

    # shellcheck disable=SC2029
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        "ubuntu@$API_IP_ADDRESS" "docker run \
       --env-file environment \
       -e DATABASE_HOST=$DATABASE_HOST \
       -e DATABASE_NAME=$DATABASE_NAME \
       -e DATABASE_USER=$DATABASE_USER \
       -e DATABASE_PASSWORD=$DATABASE_PASSWORD \
       -e ELASTICSEARCH_HOST=$ELASTICSEARCH_HOST \
       -e ELASTICSEARCH_PORT=$ELASTICSEARCH_PORT \
       -v /tmp/volumes_static:/tmp/www/static \
       --log-driver=awslogs \
       --log-opt awslogs-region=$AWS_REGION \
       --log-opt awslogs-group=data-refinery-log-group-$USER-$STAGE \
       --log-opt awslogs-stream=log-stream-api-$USER-$STAGE \
       -p 8081:8081 \
       --name=dr_api \
       -it -d $DOCKERHUB_REPO/$API_DOCKER_IMAGE /bin/sh -c /home/user/collect_and_run_uwsgi.sh"

    # shellcheck disable=SC2029
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        "ubuntu@$API_IP_ADDRESS" "docker run \
       --env-file environment \
       -e DATABASE_HOST=$DATABASE_HOST \
       -e DATABASE_NAME=$DATABASE_NAME \
       -e DATABASE_USER=$DATABASE_USER \
       -e DATABASE_PASSWORD=$DATABASE_PASSWORD \
       -e ELASTICSEARCH_HOST=$ELASTICSEARCH_HOST \
       -e ELASTICSEARCH_PORT=$ELASTICSEARCH_PORT \
       -v /tmp/volumes_static:/tmp/www/static \
       --log-driver=awslogs \
       --log-opt awslogs-region=$AWS_REGION \
       --log-opt awslogs-group=data-refinery-log-group-$USER-$STAGE \
       --log-opt awslogs-stream=log-stream-api-$USER-$STAGE \
       --name=dr_api_stats_refresh \
       -it -d $DOCKERHUB_REPO/$API_DOCKER_IMAGE python3 manage.py start_stats_refresh"

    # Don't leave secrets lying around.
    ssh -o StrictHostKeyChecking=no \
        -o ServerAliveInterval=15 \
        -o ConnectTimeout=5 \
        -i data-refinery-key.pem \
        "ubuntu@$API_IP_ADDRESS" "rm -f environment"
fi

echo "Deploy completed successfully."
