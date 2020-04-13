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
apt-get install --yes jq iotop dstat speedometer awscli docker.io chrony htop monit

ulimit -n 65536

# Find, configure and mount a free EBS volume
mkdir -p /var/ebs/

# Takes USER, STAGE
fetch_and_mount_volume () {
    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)

    # Try to mount volume 0 first, so we have one volume we know is always mounted!
    if aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=tag:Index,Values=0" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | grep 'VolumeID'; then
        EBS_VOLUME_ID=`aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=tag:Index,Values=0" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"'`
    else
        EBS_VOLUME_ID=`aws ec2 describe-volumes --filters "Name=tag:User,Values=$1" "Name=tag:Stage,Values=$2" "Name=tag:IsBig,Values=True" "Name=status,Values=available" "Name=availability-zone,Values=us-east-1a" --region us-east-1 | jq '.Volumes[0].VolumeId' | tr -d '"'`
    fi

    aws ec2 attach-volume --volume-id $EBS_VOLUME_ID --instance-id $INSTANCE_ID --device "/dev/sdf" --region ${region}
}

export STAGE=${stage}
export USER=${user}

until fetch_and_mount_volume "$USER" "$STAGE"; do
    sleep 10
done

COUNTER=0
while [  $COUNTER -lt 99 ]; do
        EBS_VOLUME_INDEX=`aws ec2 describe-volumes --filters "Name=tag:Index,Values=*" "Name=volume-id,Values=$EBS_VOLUME_ID" --query "Volumes[*].{ID:VolumeId,Tag:Tags}" --region ${region} | jq ".[0].Tag[$COUNTER].Value" | tr -d '"'`
        if echo "$EBS_VOLUME_INDEX" | egrep -q '^\-?[0-9]+$'; then
            echo "$EBS_VOLUME_INDEX is an integer!"
            break # This is a Volume Index
        else
            echo "$EBS_VOLUME_INDEX is not an integer"
        fi
        let COUNTER=COUNTER+1
done

sleep 25
# We want to mount the biggest volume that its attached to the instance
# The size of this volume can be controlled with the varialbe
# `volume_size_in_gb` from the file `variables.tf`
ATTACHED_AS=`lsblk -n --sort SIZE | tail -1 | cut -d' ' -f1`

# grep -v ext4: make sure the disk is not already formatted.
if file -s /dev/$ATTACHED_AS | grep data | grep -v ext4; then
	mkfs -t ext4 /dev/$ATTACHED_AS # This is slow
fi
mount /dev/$ATTACHED_AS /var/ebs/

chown ubuntu:ubuntu /var/ebs/
echo $EBS_VOLUME_INDEX >  /var/ebs/VOLUME_INDEX
chown ubuntu:ubuntu /var/ebs/VOLUME_INDEX

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

# Set up the Docker hung process killer
cat <<EOF >/home/ubuntu/killer.py
# Call like:
# docker ps --format 'table {{.Names}}|{{.RunningFor}}' | grep -v qn | grep -v compendia | python killer.py

import os
import sys

dockerps = sys.stdin.read()

for item in dockerps.split('\n'):
    # skip the first
    if 'NAMES' in item:
        continue
    if item == '':
        continue

    cid, time = item.split('|')
    if 'hours' not in time:
        continue

    num_hours = int(time.split(' ')[0])
    if num_hours > 1:
        print("Killing " + cid)
        os.system('docker kill ' + cid)
EOF
# Create the CW metric job in a crontab
# write out current crontab
crontab -l > tempcron
echo -e "SHELL=/bin/bash\nPATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin\n*/5 * * * * docker ps --format 'table {{.Names}}|{{.RunningFor}}' | grep -v qn | grep -v compendia | python /home/ubuntu/killer.py" >> tempcron
# install new cron file
crontab tempcron
rm tempcron

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
