#!/bin/bash -x

# This is a template for the instance-user-data.sh script for the nomad client.
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

export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get upgrade -y
apt-get install --yes jq iotop dstat speedometer awscli docker.io monit

ulimit -n 65536

# Set up the data directory
mkdir -p /var/ebs/
chown ubuntu:ubuntu /var/ebs/

# Set up the required database extensions.
# HStore allows us to treat object annotations as pseudo-NoSQL data tables.
apt-get install --yes postgresql-client-common postgresql-client
PGPASSWORD=${database_password} psql -c 'CREATE EXTENSION IF NOT EXISTS hstore;' -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}

# Change to home directory of the default user
cd /home/ubuntu

# Install, configure and launch our CloudWatch Logs agent
cat <<EOF >awslogs.conf
[general]
state_file = /var/lib/awslogs/agent-state
EOF

mkdir /var/lib/awslogs
wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
python ./awslogs-agent-setup.py --region ${region} --non-interactive --configfile awslogs.conf
echo "
/var/log/nomad_client.log {
    missingok
    notifempty
    compress
    size 20k
    daily
    maxage 3
}" >> /etc/logrotate.conf

# Docker runs out of IPv4
# Cannot specify bip option in config file because it is hardcoded in
# the startup command because docker is run by clowns.
service docker stop
nohup /usr/bin/dockerd -s overlay2 --bip=172.17.77.1/22 --log-driver=json-file --log-opt max-size=100m --log-opt max-file=3 > /var/log/docker_daemon.log &

# Output the files we need to start up Nomad and register jobs:
# (Note that the lines starting with "$" are where
#  Terraform will template in the contents of those files.)

# Create the script to install Nomad.
cat <<"EOF" > install_nomad.sh
${install_nomad_script}
EOF

# Create the Nomad Client configuration.
cat <<"EOF" > client.hcl
${nomad_client_smasher_config}
EOF

# Install Nomad
chmod +x install_nomad.sh
./install_nomad.sh

# Start the Nomad agent in client mode via Monit
echo "
#!/bin/sh
nomad status
exit \$?
" >> /home/ubuntu/nomad_status.sh
chmod +x /home/ubuntu/nomad_status.sh

echo "
#!/bin/sh
killall nomad
sleep 120
nomad agent -config /home/ubuntu/client.hcl > /var/log/nomad_client.log &
" >> /home/ubuntu/kill_restart_nomad.sh
chmod +x /home/ubuntu/kill_restart_nomad.sh
/home/ubuntu/kill_restart_nomad.sh

echo '
check program nomad with path "/bin/bash /home/ubuntu/nomad_status.sh" as uid 0 and with gid 0
    start program = "/home/ubuntu/kill_restart_nomad.sh" as uid 0 and with gid 0 with timeout 240 seconds
    if status != 0
        then restart
set daemon 300
' >> /etc/monit/monitrc

service monit restart

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
