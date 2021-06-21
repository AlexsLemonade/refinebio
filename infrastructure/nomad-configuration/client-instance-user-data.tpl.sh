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

# Configure and mount the EBS volume
mkdir -p /var/ebs/

export STAGE="${stage}"
export USER="${user}"
EBS_VOLUME_INDEX="$(wget -q -O - http://169.254.169.254/latest/meta-data/instance-id)"

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

# Kick off a script that will clean up all of our nomad jobs when this instance goes down
setsid sh clean-nomad-jobs.sh

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
