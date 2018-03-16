#!/bin/bash

# This is a template for the instance-user-data.sh script for the Nomad Server.
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


# Change to home directory of the default user
cd /home/ubuntu

# Install, configure and launch our CloudWatch Logs agent
cat <<EOF >awslogs.conf
[general]
state_file = /var/lib/awslogs/agent-state

[/var/log/nomad_server.log]
file = /var/log/nomad_server.log
log_group_name = data-refinery-log-group-${user}-${stage}
log_stream_name = log-stream-nomad-server-${server_number}-${user}-${stage}
EOF

mkdir /var/lib/awslogs
wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
python ./awslogs-agent-setup.py --region ${region} --non-interactive --configfile awslogs.conf
# Rotate the logs, delete after 3 days.
echo "
/var/log/nomad_server.log {
    missingok
    notifempty
    compress
    size 20k
    daily
    maxage 3
}" >> /etc/logrotate.conf

# Output the files we need to start up Nomad.
# Note that the lines starting with "$" are where
# Terraform will template in the contents of those files.

# Create the script to install Nomad.
cat <<"EOF" > install_nomad.sh
${install_nomad_script}
EOF

# Create the Nomad Server configuration.
cat <<"EOF" > server.hcl
${nomad_server_config}
EOF

# Install Nomad.
chmod +x install_nomad.sh
./install_nomad.sh

# Start the Nomad agent in server mode.
nohup nomad agent -config server.hcl > /var/log/nomad_server.log &

# Give the Nomad server time to start up.
sleep 30

# Determine the IP address of this machine so we know where the Nomad
# server's HTTP API is addressed.
IP_ADDRESS=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

# Register job specs.
nomad_job_specs=nomad-job-specs/*
for nomad_job_spec in $nomad_job_specs; do
    echo "registering $nomad_job_spec"
    nomad run -address http://$IP_ADDRESS:4646 $nomad_job_spec
done

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi