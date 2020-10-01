#!/bin/bash -x

# This is a template for the instance-user-data.sh script for the nomad client.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from the
# project into terraform variables, and then template them into the following
# script. These will then be written out to files so that they can be used
# locally. This means that any variable referenced using `{}` is NOT a shell
# variable, it is a template variable for Terraform to fill in. DO NOT treat
# them as normal shell variables.

# A more ideal solution than this would be if we could just give AWS a list of
# files to put onto the instance for us, but they only give us this one script
# to do it with. Nomad has file provisioners which will put files onto the
# instance after it starts up, but those run after this script runs.

export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get upgrade -y
apt-get install --yes jq iotop dstat speedometer awscli docker.io

ulimit -n 65536

# Find, configure and mount the EBS volume
mkdir -p /var/ebs/

# Takes USER, STAGE
fetch_and_mount_volume () {
    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)

    EBS_VOLUME_ID=$(aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"')

    aws ec2 attach-volume --volume-id "$EBS_VOLUME_ID" --instance-id "$INSTANCE_ID" --device "/dev/sdf" --region "${region}"
}

export STAGE="${stage}"
export USER="${user}"
EBS_VOLUME_INDEX="$(wget -q -O - http://169.254.169.254/latest/meta-data/instance-id)"

until fetch_and_mount_volume "$USER" "$STAGE"; do
    sleep 10
done

# We want to mount the biggest volume that its attached to the instance
# The size of this volume can be controlled with the varialbe
# `smasher_volume_size_in_gb` from the file `variables.tf`
ATTACHED_AS=$(lsblk -n --sort SIZE | tail -1 | cut -d' ' -f1)

# # grep -v ext4: make sure the disk is not already formatted.
if file -s "/dev/$ATTACHED_AS" | grep data | grep -v ext4; then
	mkfs -t ext4 "/dev/$ATTACHED_AS" # This is slow
fi
mount "/dev/$ATTACHED_AS" /var/ebs/

chown ubuntu:ubuntu /var/ebs/
echo "$EBS_VOLUME_INDEX" > /var/ebs/VOLUME_INDEX
chown ubuntu:ubuntu /var/ebs/VOLUME_INDEX


# Set up the required database extensions.
# HStore allows us to treat object annotations as pseudo-NoSQL data tables.
apt-get install --yes postgresql-client-common postgresql-client
PGPASSWORD=${database_password} psql -c 'CREATE EXTENSION IF NOT EXISTS hstore;' -h "${database_host}" -p 5432 -U "${database_user}" -d "${database_name}"

# Change to home directory of the default user
cd /home/ubuntu || exit

# Install, configure and launch our CloudWatch Logs agent
cat <<EOF >awslogs.conf
[general]
state_file = /var/lib/awslogs/agent-state
EOF

mkdir /var/lib/awslogs
wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
python3.5 ./awslogs-agent-setup.py --region "${region}" --non-interactive --configfile awslogs.conf
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

# Start the Nomad agent in client mode via systemd
cat <<"EOF" > /etc/systemd/system/nomad-client.service
${nomad_client_service}
EOF

systemctl enable nomad-client.service
systemctl start nomad-client.service

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
