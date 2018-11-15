#!/bin/bash -e

print_description() {
    echo 'This script can be used to deploy and update a `refine.bio` instance stack.'
    echo 'It will create all of the AWS infrasctructure (roles/instances/db/network/etc),'
    echo 'open an ingress, kill all running Nomad jobs, perform a database migration,'
    echo 're-define and re-register Nomad job specifications, and finally close the'
    echo 'ingress. This can be run from a CI/CD machine or a local dev box.'
    echo 'This script must be run from /infrastructure!'
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
    echo '   If unset, defaults to the value in `infrastructure/environments/$env`, which is "ccdlstaging"'
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

if [[ -z $TF_VAR_region ]]; then
    TF_VAR_region=us-east-1
fi

# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  $NOMAD_ADDR/v1/status/leader)
}

# We have terraform output environment variables via a single output
# variable, which we then read in as json using the command line tool
# `jq`, so that we can use them via bash.
format_environment_variables () {
  json_env_vars=$(terraform output -json environment_variables | jq -c '.value[]')
  for row in $json_env_vars; do
      env_var_assignment=$(echo $row | jq -r ".name")=$(echo $row | jq -r ".value")
      export $env_var_assignment
      echo $env_var_assignment >> prod_env
  done
}

# Load $ALL_CCDL_IMAGES and helper functions
source ../common.sh
# Make our IP address known to terraform.
export TF_VAR_host_ip=`dig +short myip.opendns.com @resolver1.opendns.com`

for IMAGE in $ALL_CCDL_IMAGES; do
    # For each image we need to set the env var that is used by our
    # scripts and the env var that gets picked up by terraform because
    # it is preceeded with TF_VAR.
    IMAGE_UPPER=$IMAGE | tr a-z A-Z
    export ${IMAGE_UPPER}_DOCKER_IMAGE=dr_$IMAGE:$SYSTEM_VERSION
    export TF_VAR_${IMAGE}_docker_image=dr_$IMAGE:$SYSTEM_VERSION
done

# Copy ingress config to top level so it can be applied.
cp deploy/ci_ingress.tf .

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

    if [[ ! -z $CIRCLE_BUILD_NUM ]]; then
        # Make sure we can't expose secrets in circleci
        terraform apply -var-file=environments/$env.tfvars -auto-approve > /dev/null
    else
        terraform apply -var-file=environments/$env.tfvars -auto-approve
    fi
fi

# We have to do this once before the initial deploy..
rm -f prod_env
format_environment_variables

../format_nomad_with_env.sh -p api -e $env -o $(pwd)/api-configuration/
../format_nomad_with_env.sh -p foreman -e $env -o $(pwd)/foreman-configuration/

if [[ -z $ran_init_build ]]; then
    # Open up ingress to AWS for Circle, stop jobs, migrate DB.
    echo "Deploying with ingress.."

    # Output the plan for debugging deployments later.
    # Until terraform plan supports -var-file the plan is wrong.
    # terraform plan

    if [[ ! -z $CIRCLE_BUILD_NUM ]]; then
        # Make sure we can't expose secrets in circleci
        terraform apply -var-file=environments/$env.tfvars -auto-approve > /dev/null
    else
        terraform apply -var-file=environments/$env.tfvars -auto-approve
    fi
fi

# Find address of Nomad server.
export NOMAD_LEAD_SERVER_IP=`terraform output nomad_server_1_ip`

export NOMAD_ADDR=http://$NOMAD_LEAD_SERVER_IP:4646

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
echo "Killing base jobs.. (this takes a while..)"
if [[ $(nomad status) != "No running jobs" ]]; then
    for job in $(nomad status | grep running | awk {'print $1'} || grep --invert-match /)
    do
        # '|| true' so that if a job is garbage collected before we can remove it the error
        # doesn't interrupt our deploy.
        nomad stop -purge -detach $job > /dev/null || true &
    done
fi

# Wait to make sure that all base jobs are killed so no new jobs can
# be queued while we kill the parameterized Nomad jobs.
wait $(jobs -p)

# Kill parameterized Nomad Jobs so no jobs will be running when we
# apply migrations.
echo "Killing dispatch jobs.. (this also takes a while..)"
if [[ $(nomad status) != "No running jobs" ]]; then
    counter=0
    for job in $(nomad status | awk {'print $1'} || grep /)
    do
        # Skip the header row for jobs.
        if [ $job != "ID" ]; then
            # '|| true' so that if a job is garbage collected before we can remove it the error
            # doesn't interrupt our deploy.
            nomad stop -purge -detach $job > /dev/null || true &
            counter=$((counter+1))
        fi

        # Wait for all the jobs to stop every 100 so we don't knock
        # over the deploy box if there are 1000's.
        if [[ "$counter" -gt 100 ]]; then
            wait $(jobs -p)
            counter=0
        fi
    done
fi

# Wait for any remaining jobs to all die.
wait $(jobs -p)

# Make sure that prod_env is empty since we are only appending to it.
# prod_env is a temporary file we use to pass environment variables to
# `docker run` commands when running migrations.
rm -f prod_env

# (cont'd) ..and once again after the update when this is re-run.
format_environment_variables

# Get an image to run the migrations with.
docker pull $DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE

# Migrate auth.
# Override database settings so we connect directly to the RDS instance.
docker run \
       --env-file prod_env \
       --env DATABASE_HOST=$RDS_HOST \
       --env DATABASE_PORT=$DATABASE_HIDDEN_PORT \
       --env RUNNING_IN_CLOUD=False \
       $DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE python3 manage.py migrate auth

# Apply general migrations.
docker run \
       --env-file prod_env \
       --env DATABASE_HOST=$RDS_HOST \
       --env DATABASE_PORT=$DATABASE_HIDDEN_PORT \
       --env RUNNING_IN_CLOUD=False \
       $DOCKERHUB_REPO/$FOREMAN_DOCKER_IMAGE python3 manage.py migrate

# Template the environment variables for production into the Nomad Job
# specs and API confs.
mkdir -p nomad-job-specs
../format_nomad_with_env.sh -p workers -e $env -o $(pwd)/nomad-job-specs
../format_nomad_with_env.sh -p surveyor -e $env -o $(pwd)/nomad-job-specs

# API and foreman aren't run as nomad jobs, but the templater still works.
../format_nomad_with_env.sh -p foreman -e $env -o $(pwd)/foreman-configuration
../format_nomad_with_env.sh -p api -e $env -o $(pwd)/api-configuration/


# Don't leave secrets lying around!
rm -f prod_env

# Re-register Nomad jobs.
echo "Registering new job specifications.."
nomad_job_specs=nomad-job-specs/*
for nomad_job_spec in $nomad_job_specs; do
    nomad run $nomad_job_spec &
done
echo "Job registrations have been fired off."

# Ensure the latest image version is being used for the Foreman
terraform taint aws_instance.foreman_server_1

# Remove the ingress config so the next `terraform apply` will remove
# access for Circle.
echo "Removing ingress.."
rm ci_ingress.tf

if [[ ! -z $CIRCLE_BUILD_NUM ]]; then
    # Make sure we can't expose secrets in circleci
    terraform apply -var-file=environments/$env.tfvars -auto-approve > /dev/null
else
    terraform apply -var-file=environments/$env.tfvars -auto-approve
fi

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
                        -i data-refinery-key.pem \
                        ubuntu@$API_IP_ADDRESS  "docker ps -a" | grep dr_api || echo "")

# If $container_running is empty, then it's because the container isn't running.
# If the container isn't running, then it's because the instance is spinning up.
# The container will be started by the API's init script, so no need to do anything more.
# However if $container_running isn't empty then we need to stop and restart it.
if [[ ! -z $container_running ]]; then
    echo "Restarting API with latest image."

    ssh -o StrictHostKeyChecking=no \
        -i data-refinery-key.pem \
        ubuntu@$API_IP_ADDRESS  "docker pull $DOCKERHUB_REPO/$API_DOCKER_IMAGE"

    ssh -o StrictHostKeyChecking=no \
        -i data-refinery-key.pem \
        ubuntu@$API_IP_ADDRESS "docker rm -f dr_api"

    scp -o StrictHostKeyChecking=no \
        -i data-refinery-key.pem \
        api-configuration/environment ubuntu@$API_IP_ADDRESS:/home/ubuntu/environment

    ssh -o StrictHostKeyChecking=no \
        -i data-refinery-key.pem \
        ubuntu@$API_IP_ADDRESS "docker run \
       --env-file environment \
       -e DATABASE_HOST=$DATABASE_HOST \
       -e DATABASE_NAME=$DATABASE_NAME \
       -e DATABASE_USER=$DATABASE_USER \
       -e DATABASE_PASSWORD=$DATABASE_PASSWORD \
       -v /tmp/volumes_static:/tmp/www/static \
       --log-driver=awslogs \
       --log-opt awslogs-region=$REGION \
       --log-opt awslogs-group=data-refinery-log-group-$USER-$STAGE \
       --log-opt awslogs-stream=log-stream-api-$USER-$STAGE \
       -p 8081:8081 \
       --name=dr_api \
       -it -d $DOCKERHUB_REPO/$API_DOCKER_IMAGE /bin/sh -c /home/user/collect_and_run_uwsgi.sh"

    # Don't leave secrets lying around.
    ssh -o StrictHostKeyChecking=no \
        -i data-refinery-key.pem \
        ubuntu@$API_IP_ADDRESS "rm -f environment"
fi

echo "Deploy completed successfully."
