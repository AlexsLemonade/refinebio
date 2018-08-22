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

# Mount the EFS.
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get upgrade -y
apt-get install --yes nfs-common jq iotop dstat

# EFS
# mkdir -p /var/efs/
# echo "${file_system_id}.efs.${region}.amazonaws.com:/ /var/efs/ nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2 0 0" >> /etc/fstab
# mount -a -t nfs4
# chown ubuntu:ubuntu /var/efs/

# Find, configure and mount a free EBS volume
mkdir -p /var/efs/
EBS_VOLUME_ID=`aws ec2 describe-volumes --filters "Name=tag:User,Values=${user}" "Name=tag:Stage,Values=${stage}" "Name=tag:IsBig,Values=True" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"'`

COUNTER=0
while [  $COUNTER -lt 5 ]; do
        EBS_VOLUME_INDEX=`aws ec2 describe-volumes --filters "Name=tag:Index,Values=*" "Name=volume-id,Values=$EBS_VOLUME_ID" --query "Volumes[*].{ID:VolumeId,Tag:Tags}" --region us-east-1 | jq ".[0].Tag[$COUNTER].Value" | tr -d '"'`
        if echo "$EBS_VOLUME_INDEX" | egrep -q '^\-?[0-9]+$'; then
            echo "$EBS_VOLUME_INDEX is an integer!"
            break # This is a Volume Index
        else
            echo "$EBS_VOLUME_INDEX is not an integer"
        fi
        let COUNTER=COUNTER+1
done

INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)
aws ec2 attach-volume --volume-id $EBS_VOLUME_ID --instance-id $INSTANCE_ID --device "/dev/sdf" --region us-east-1 # <_<
sleep 15
ATTACHED_AS=`lsblk -n | grep T | cut -d' ' -f1`
FILE_RESULT=`file -s /dev/$ATTACHED_AS`

if file -s /dev/$ATTACHED_AS | grep data; then
	mkfs -t ext4 /dev/$ATTACHED_AS # This is slow
fi
mount /dev/$ATTACHED_AS /var/efs/
chown ubuntu:ubuntu /var/efs/
echo $EBS_VOLUME_INDEX >  /var/efs/VOLUME_INDEX
chown ubuntu:ubuntu /var/efs/VOLUME_INDEX

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

[/var/log/nomad_client.log]
file = /var/log/nomad_client.log
log_group_name = data-refinery-log-group-${user}-${stage}
log_stream_name = log-stream-nomad-client-${user}-${stage}
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

# Output the files we need to start up Nomad and register jobs:
# (Note that the lines starting with "$" are where
#  Terraform will template in the contents of those files.)

# Create the script to install Nomad.
cat <<"EOF" > install_nomad.sh
${install_nomad_script}
EOF

# Create the Nomad Client configuration.
cat <<"EOF" > client.hcl
${nomad_client_config}
EOF
# Make the client.meta.volume_id is set to waht we just mounted
sed -i "s/REPLACE_ME/$EBS_VOLUME_INDEX/" client.hcl

# Create a directory for docker to use as a volume.
mkdir /home/ubuntu/docker_volume
chmod a+rwx /home/ubuntu/docker_volume

# Install Nomad
chmod +x install_nomad.sh
./install_nomad.sh

# Start the Nomad agent in client mode.
nomad agent -config client.hcl > /var/log/nomad_client.log &

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
