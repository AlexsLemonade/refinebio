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

# Create and install SSL Certificate
apt-get install software-properties-common
add-apt-repository ppa:certbot/certbot
apt-get update
apt-get install python-certbot-nginx

certbot --nginx -d api.staging.refine.bio -n --agree-tos -m g3w4k4t5n3s7p7v8@alexslemonade.slack.com

# Let's Encrypt certs are eligible for renewal 30 days before they
# expire so try to renew every 29 days and we'll be sure to hit the
# window. (They are valid for 90 days and renewing too early has no
# consequences other than potentially hitting a rate limit if you do
# it too often, which this won't.)
# See https://certbot.eff.org/docs/using.html#renewing-certificates
# The command to create this crontab was made with help from:
# https://stackoverflow.com/a/610860/6095378
# and
# https://stackoverflow.com/a/550808/6095378
echo '0 0 */29 * * root certbot renew' | sudo tee --append /etc/cron.d/renew_ssl_cert

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

STATIC_VOLUMES=/tmp/volumes_static
mkdir -p /tmp/volumes_static
chmod a+rwx /tmp/volumes_static

# These database values are created after TF
# is run, so we have to pass them in programatically
docker run \
       --env-file environment \
       -e DATABASE_HOST=${database_host} \
       -e DATABASE_NAME=${database_name} \
       -e DATABASE_USER=${database_user} \
       -e DATABASE_PASSWORD=${database_password} \
       -v "$STATIC_VOLUMES":/tmp/www/static \
       --log-driver=awslogs \
       --log-opt awslogs-region=${region} \
       --log-opt awslogs-group=${log_group} \
       --log-opt awslogs-stream=${log_stream} \
       -p 8081:8081 \
       -it -d ${dockerhub_repo}/${api_docker_image} /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"


# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
