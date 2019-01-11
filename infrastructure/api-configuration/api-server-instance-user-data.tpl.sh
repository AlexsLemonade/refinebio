#!/bin/bash

# This is a template for the instance-user-data.sh script for the API Server.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from
# the project into terraform variables, and then template them into
# the following script. These will then be written out to files
# so that they can be used locally.

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

if [[ ${stage} == "staging" || ${stage} == "prod" ]]; then
    # Create and install SSL Certificate for the API.
    # Only necessary on staging and prod.
    # We cannot use ACM for this because *.bio is not a Top Level Domain that Route53 supports.
    apt-get install -y software-properties-common
    add-apt-repository ppa:certbot/certbot
    apt-get update
    apt-get install -y python-certbot-nginx

    # g3w4k4t5n3s7p7v8@alexslemonade.slack.com is the email address we
    # have configured to forward mail to the #teamcontact channel in
    # slack. Certbot will use it for "important account
    # notifications".
    # In the future, if we ever hit the 5-deploy-a-week limit, changing one of these lines to:
    # certbot --nginx -d api.staging.refine.bio -d api2.staging.refine.bio -n --agree-tos --redirect -m g3w4k4t5n3s7p7v8@alexslemonade.slack.com
    # will circumvent certbot's limit because the 5-a-week limit only applies to the
    # "same set of domains", so by changing that set we get to use the 20-a-week limit.
    if [[ ${stage} == "staging" ]]; then
        # The certbot challenge cannot be completed until the aws_lb_target_group_attachment resources are created.
        sleep 180
        certbot --nginx -d api.staging.refine.bio -n --agree-tos --redirect -m g3w4k4t5n3s7p7v8@alexslemonade.slack.com
    elif [[ ${stage} == "prod" ]]; then
        # The certbot challenge cannot be completed until the aws_lb_target_group_attachment resources are created.
        sleep 180
        certbot --nginx -d api.refine.bio -n --agree-tos --redirect -m g3w4k4t5n3s7p7v8@alexslemonade.slack.com
    fi
fi

# Install, configure and launch our CloudWatch Logs agent
cat <<EOF >awslogs.conf
[general]
state_file = /var/lib/awslogs/agent-state

[/tmp/error.log]
file = /tmp/error.log
log_group_name = ${log_group}
log_stream_name = log-stream-api-nginx-error-${user}-${stage}

[/tmp/access.log]
file = /tmp/access.log
log_group_name = ${log_group}
log_stream_name = log-stream-api-nginx-access-${user}-${stage}

EOF

mkdir /var/lib/awslogs
wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
python ./awslogs-agent-setup.py --region ${region} --non-interactive --configfile awslogs.conf
# Rotate the logs, delete after 3 days.
echo "
/tmp/error.log {
    missingok
    notifempty
    compress
    size 20k
    daily
    maxage 3
}" >> /etc/logrotate.conf
echo "
/tmp/access.log {
    missingok
    notifempty
    compress
    size 20k
    daily
    maxage 3
}" >> /etc/logrotate.conf

# Install our environment variables
cat <<"EOF" > environment
${api_environment}
EOF

chown -R ubuntu /home/ubuntu

STATIC_VOLUMES=/tmp/volumes_static
mkdir -p /tmp/volumes_static
chmod a+rwx /tmp/volumes_static

# Pull the API image.
docker pull ${dockerhub_repo}/${api_docker_image}

# These database values are created after TF
# is run, so we have to pass them in programatically
docker run \
       --env-file environment \
       -e DATABASE_HOST=${database_host} \
       -e DATABASE_NAME=${database_name} \
       -e DATABASE_USER=${database_user} \
       -e DATABASE_PASSWORD=${database_password} \
       -e ELASTICSEARCH_HOST=${elasticsearch_host} \
       -e ELASTICSEARCH_PORT=${elasticsearch_port} \
       -v "$STATIC_VOLUMES":/tmp/www/static \
       --log-driver=awslogs \
       --log-opt awslogs-region=${region} \
       --log-opt awslogs-group=${log_group} \
       --log-opt awslogs-stream=${log_stream} \
       -p 8081:8081 \
       --name=dr_api \
       -it -d ${dockerhub_repo}/${api_docker_image} /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"

# Let's index now too
docker exec -it dr_api python3 manage.py search_index --rebuild -f

# Don't leave secrets lying around.
rm -f environment

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
