#!/bin/bash

# This is a template for the instance-user-data.sh script for the Foreman Server.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from
# the project into terraform variables, and then template them into
# the following script. These will then be written out to files
# so that they can be used locally.

# Change to home directory of the default user
cd /home/ubuntu

# Install our environment variables
cat <<"EOF" > environment
${foreman_environment}
EOF

# These database values are created after TF
# is run, so we have to pass them in programatically
docker run \
       --env-file environment \
       -e DATABASE_HOST=${database_host} \
       -e DATABASE_NAME=${database_name} \
       -e DATABASE_USER=${database_user} \
       -e DATABASE_PASSWORD=${database_password} \
       --add-host=nomad:${nomad_lead_server_ip}\
       --log-driver=awslogs \
       --log-opt awslogs-region=${region} \
       --log-opt awslogs-group=${log_group} \
       --log-opt awslogs-stream=log-stream-foreman-${user}-${stage} \
       -it -d ${dockerhub_repo}/${foreman_docker_image} python3 manage.py retry_jobs


# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
