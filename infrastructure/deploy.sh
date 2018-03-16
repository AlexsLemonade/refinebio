#!/bin/bash -e

# This script can be used to deploy and update a `refine.bio` instance stack.
# It will create all of the AWS infrasctructure (roles/instances/db/network/etc),
# open an ingress, kill all running Nomad jobs, perform a database migration,
# re-define and re-register Nomad job specifications, and finally close the
# ingress. This can be run from a CI/CD machine or a local dev box.
# This script must be run from /infrastructure!

# Make our IP address known to terraform.
source ../common.sh
export TF_VAR_host_ip=`dig +short myip.opendns.com @resolver1.opendns.com`

# Copy ingress config to top level so it can be applied.
cp deploy/ci_ingress.tf .

# Output the plan for debugging deployments later.
terraform plan

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
        nomad stop -purge -detach $job > /dev/null
    done
fi

# Kill parameterized Nomad Jobs so no jobs will be running when we
# apply migrations.
echo "Killing dispatch jobs... (This may take a while.)"
if [[ $(nomad status) != "No running jobs" ]]; then
    for job in $(nomad status | awk {'print $1'} || grep /)
    do
        # Skip the header row for jobs.
        if [ $job != "ID" ]; then
            nomad stop -purge -detach $job > /dev/null
        fi
    done
fi

# Make sure that prod_env is empty since we are only appending to it.
# prod_env is a temporary file we use to pass environment variables to
# `docker run` commands when running migrations.
rm -f prod_env

# We have terraform output environment variables via a single output
# variable, which we then read in as json using the command line tool
# `jq`, so that we can use them via bash.
for row in $(terraform output -json environment_variables | jq -c '.value[]'); do
    env_var_assignment=$(echo $row | jq -r ".name")=$(echo $row | jq -r ".value")
    export $env_var_assignment
    echo $env_var_assignment >> prod_env
done

# Get an image to run the migrations with.
docker pull $FOREMAN_DOCKER_IMAGE

# Migrate auth.
docker run \
       --volume migrations \
       --env-file prod_env \
       $FOREMAN_DOCKER_IMAGE migrate auth

# Apply general migrations.
docker run \
       --volume migrations \
       --env-file prod_env \
       $FOREMAN_DOCKER_IMAGE migrate

# Don't leave secrets lying around!
rm prod_env

# Template the environment variables for production into the Nomad Job
# specs.
mkdir -p nomad-job-specs
../format_nomad_with_env.sh -p workers -e prod -o $(pwd)/nomad-job-specs/
../format_nomad_with_env.sh -p foreman -e prod -o $(pwd)/nomad-job-specs/

# Re-register Nomad jobs.
echo "Registering new job specifications.."
nomad_job_specs=nomad-job-specs/*
for nomad_job_spec in $nomad_job_specs; do
    echo "Registering $nomad_job_spec"
    nomad run $nomad_job_spec
done

# Remove the ingress config so the next `terraform apply` will remove
# access for Circle.
echo "Removing ingress.."
rm ci_ingress.tf
terraform apply -auto-approve
