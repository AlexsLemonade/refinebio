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
echo "
#!/bin/sh
docker rm -f \$(docker ps --quiet --all) || true
docker run \\
       --env-file /home/ubuntu/environment \\
       -e DATABASE_HOST=${database_host} \\
       -e DATABASE_NAME=${database_name} \\
       -e DATABASE_USER=${database_user} \\
       -e DATABASE_PASSWORD=${database_password} \\
       -e ELASTICSEARCH_HOST=${elasticsearch_host} \\
       -e ELASTICSEARCH_PORT=${elasticsearch_port} \\
       -v /tmp:/tmp \\
       --add-host=nomad:${nomad_lead_server_ip} \\
       --log-driver=awslogs \\
       --log-opt awslogs-region=${region} \\
       --log-opt awslogs-group=${log_group} \\
       --log-opt awslogs-stream=log-stream-foreman-${user}-${stage} \\
       --name=dr_foreman \\
       -it -d ${dockerhub_repo}/${foreman_docker_image} python3 manage.py retry_jobs
" >> /home/ubuntu/run_foreman.sh
chmod +x /home/ubuntu/run_foreman.sh
/home/ubuntu/run_foreman.sh

# The foreman instance is used for running various management
# commands. This script provides an easy way to do that.
echo "
#!/bin/sh

# This script should be used by passing a management command name as
# the first argument followed by any additional arguments to that
# command.

docker run \\
       --env-file /home/ubuntu/environment \\
       -e DATABASE_HOST=${database_host} \\
       -e DATABASE_NAME=${database_name} \\
       -e DATABASE_USER=${database_user} \\
       -e DATABASE_PASSWORD=${database_password} \\
       -e ELASTICSEARCH_HOST=${elasticsearch_host} \\
       -e ELASTICSEARCH_PORT=${elasticsearch_port} \\
       -v /tmp:/tmp \\
       --add-host=nomad:${nomad_lead_server_ip} \\
       -it -d ${dockerhub_repo}/${foreman_docker_image} python3 manage.py \$@
" >> /home/ubuntu/run_management_command.sh
chmod +x /home/ubuntu/run_management_command.sh

# Start the Nomad agent in server mode via Monit
apt-get -y update
apt-get -y install monit htop

date +%s > /tmp/foreman_last_time
chown ubuntu:ubuntu /tmp/foreman_last_time
echo '
#!/bin/sh
lasttime=$(</tmp/foreman_last_time);
nowtime=`date +%s`;
((difftime = $nowtime - $lasttime));
if (( $difftime > 1800 )); then
  exit 1;
fi
exit 0;
' >> /home/ubuntu/foreman_status.sh
chmod +x /home/ubuntu/foreman_status.sh

echo '
check program foreman with path "/bin/bash /home/ubuntu/foreman_status.sh" as uid 0 and with gid 0
    start program = "/bin/bash /home/ubuntu/run_foreman.sh" as uid 0 and with gid 0
    if status != 0
        then restart
set daemon 900
' >> /etc/monit/monitrc

service monit restart

docker run \
       --env-file /home/ubuntu/environment \
       -e DATABASE_HOST=${database_host} \
       -e DATABASE_NAME=${database_name} \
       -e DATABASE_USER=${database_user} \
       -e DATABASE_PASSWORD=${database_password} \
       -v /tmp:/tmp \
       --add-host=nomad:${nomad_lead_server_ip} \
       --log-driver=awslogs \
       --log-opt awslogs-region=${region} \
       --log-opt awslogs-group=${log_group} \
       --log-opt awslogs-stream=log-stream-foreman-${user}-${stage} \
       --name=job_filler \
       -it -d ${dockerhub_repo}/${foreman_docker_image} python3 manage.py create_missing_processor_jobs

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
