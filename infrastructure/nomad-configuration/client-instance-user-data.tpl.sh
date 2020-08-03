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
apt-get install --yes jq iotop dstat speedometer awscli docker.io chrony htop

ulimit -n 65536

# Find, configure and mount a free EBS volume
mkdir -p /var/ebs/

# Takes USER, STAGE
# fetch_and_mount_volume () {
#     INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)

#     # Try to mount volume 0 first, so we have one volume we know is always mounted!
#     if aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=tag:Index,Values=0" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | grep 'VolumeID'; then
#         EBS_VOLUME_ID=`aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=tag:Index,Values=0" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"'`
#     else
#         EBS_VOLUME_ID=`aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"'`
#     fi

#     aws ec2 attach-volume --volume-id $EBS_VOLUME_ID --instance-id $INSTANCE_ID --device "/dev/sdf" --region ${region}
# }

export STAGE="${stage}"
export USER="${user}"
EBS_VOLUME_INDEX="$(wget -q -O - http://169.254.169.254/latest/meta-data/instance-id)"

# until fetch_and_mount_volume "$USER" "$STAGE"; do
#     sleep 10
# done

# COUNTER=0
# while [  $COUNTER -lt 99 ]; do
#         EBS_VOLUME_INDEX=`aws ec2 describe-volumes --filters "Name=tag:Index,Values=*" "Name=volume-id,Values=$EBS_VOLUME_ID" --query "Volumes[*].{ID:VolumeId,Tag:Tags}" --region ${region} | jq ".[0].Tag[$COUNTER].Value" | tr -d '"'`
#         if echo "$EBS_VOLUME_INDEX" | egrep -q '^\-?[0-9]+$'; then
#             echo "$EBS_VOLUME_INDEX is an integer!"
#             break # This is a Volume Index
#         else
#             echo "$EBS_VOLUME_INDEX is not an integer"
#         fi
#         let COUNTER=COUNTER+1
# done

# # Only a single volume for now
# export EBS_VOLUME_INDEX=0

# sleep 25
# # We want to mount the biggest volume that its attached to the instance
# # The size of this volume can be controlled with the varialbe
# # `volume_size_in_gb` from the file `variables.tf`
# ATTACHED_AS=`lsblk -n --sort SIZE | tail -1 | cut -d' ' -f1`

# # grep -v ext4: make sure the disk is not already formatted.
# if file -s /dev/$ATTACHED_AS | grep data | grep -v ext4; then
# 	mkfs -t ext4 /dev/$ATTACHED_AS # This is slow
# fi
# mount /dev/$ATTACHED_AS /var/ebs/

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
python ./awslogs-agent-setup.py --region "${region}" --non-interactive --configfile awslogs.conf
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

# Create the Nomad Client configuration.
cat <<"EOF" > client.hcl
${nomad_client_config}
EOF
# Make sure the client.meta.volume_id is set to what we just mounted
sed -i "s/REPLACE_ME/$EBS_VOLUME_INDEX/" client.hcl

# Create a directory for docker to use as a volume.
mkdir /home/ubuntu/docker_volume
chmod a+rwx /home/ubuntu/docker_volume

# Output the files we need to start up Nomad and register jobs. These are
# stored gzipped at the end of this script
sed '1,/^#__EOF__$/d' "$0" | base64 -d | tar xzv

# Start the Nomad agent in client mode via systemd
sudo mv nomad-client.service /etc/systemd/system/
systemctl enable nomad-client.service
systemctl start nomad-client.service

# Sleep for a bit so nomad can start
sleep 30

# Start the nomad jobs that are templated for this instance
for nomad_job_spec in nomad-job-specs/*; do
    sed -i 's/__REPLACE_ME__/'"$EBS_VOLUME_INDEX"'/g' "$nomad_job_spec"
    nomad run "$nomad_job_spec"
done
rm -r nomad-job-specs

# Set up the Docker hung process killer
# cat <<EOF >/home/ubuntu/killer.py
# # Call like:
# # docker ps --format 'table {{.Names}}|{{.RunningFor}}' | grep -v qn | grep -v compendia | python killer.py

# import os
# import sys

# dockerps = sys.stdin.read()

# for item in dockerps.split('\n'):
#     # skip the first
#     if 'NAMES' in item:
#         continue
#     if item == '':
#         continue

#     cid, time = item.split('|')
#     if 'hours' not in time:
#         continue

#     num_hours = int(time.split(' ')[0])
#     if num_hours > 1:
#         print("Killing " + cid)
#         os.system('docker kill ' + cid)
# EOF
# # Create the CW metric job in a crontab
# # write out current crontab
# crontab -l > tempcron
# echo -e "SHELL=/bin/bash\nPATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin\n*/5 * * * * docker ps --format 'table {{.Names}}|{{.RunningFor}}' | grep -v qn | grep -v compendia | python /home/ubuntu/killer.py" >> tempcron
# # install new cron file
# crontab tempcron
# rm tempcron

# Set up the AWS NTP
# via https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/set-time.html#configure_ntp
echo 'server 169.254.169.123 prefer iburst' | cat - /etc/chrony/chrony.conf > temp && mv temp /etc/chrony/chrony.conf
/etc/init.d/chrony restart

# Delete the cloudinit and syslog in production.
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
