#!/bin/bash

# This is a template for the instance-user-data.sh script for the Lead Nomad Server.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Nomad, which will read files from
# the project into terraform variables, and then template them into
# the following Here Documents. These will then be written out to files
# so that they can be used. A more ideal solution than this would be if
# we could just give AWS a list of files to put onto the instance for us,
# but they only give us this one script to do it with. Nomad has file
# provisioners which will put files onto the instance after it starts up,
# but those run after this script runs.

# This template varies from nomad-server-instance-user-data.tpl.sh in
# only a few ways.  Because this will be run first, it does not need
# to know the IP of any other server. Instead its IP will be used by
# the other Nomad Servers to join the Raft. Additionally since we only
# need to register the Nomad Jobs once so we do it in this script
# since it will only be run once in total (the other server startup
# script will be run twice, once by each of the other Nomad Servers).


# Change to home directory of the default user
cd /home/ubuntu

# Install and configure Nginx.
cat <<"EOF" > nginx.conf
${nginx_config}
EOF
apt-get update -y
apt-get install nginx -y
cp nginx.conf /etc/nginx/nginx.conf
service nginx restart

# Install, configure and launch our CloudWatch Logs agent
# cat <<EOF >awslogs.conf
# [general]
# state_file = /var/lib/awslogs/agent-state

# [/var/log/docker_api.log]
# file = /var/log/nomad_server.log
# log_group_name = data-refinery-log-group-${user}-${stage}
# log_stream_name = log-stream-api-server-${user}-${stage}
# EOF

# mkdir /var/lib/awslogs
# wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
# python ./awslogs-agent-setup.py --region ${region} --non-interactive --configfile awslogs.conf
# # Rotate the logs, delete after 3 days.
# echo "
# /var/log/docker_api.log {
#     missingok
#     notifempty
#     compress
#     size 20k
#     daily
#     maxage 3
# }" >> /etc/logrotate.conf

# Install our environment variables
cat <<"EOF" > environment
${api_environment}
EOF

STATIC_VOLUMES=/tmp/volumes_static
mkdir -p /tmp/volumes_static
chmod a+rwx /tmp/volumes_static
docker run \
       --env-file environment \
       -e DATABASE_HOST=${database_host} \
       -e DATABASE_NAME=${database_name} \
       -e DATABASE_USER=${database_user} \
       -e DATABASE_PASSWORD=${database_password} \
       -v "$STATIC_VOLUMES":/tmp/www/static \
       -p 8081:8081 \
       -it -d ${api_docker_image} /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"


# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi

# docker logs -f dr_api_prod2 > /var/log/docker_api.log &