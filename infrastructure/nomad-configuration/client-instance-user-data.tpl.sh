#!/bin/bash -ex

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
apt-get install --yes nfs-common
mkdir -p /var/efs/
chown ubuntu:ubuntu /var/efs/
echo "${file_system_id}.efs.${region}.amazonaws.com:/ /var/efs/ nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2 0 0" >> /etc/fstab
mount -a -t nfs4

# Set up the database
# XXX - Failure here is catastrphic - even normally "safe" operations like apply permissions to a user that already has them
# XXX - We need to have way to do this "manually" for staging and automatically for prod, all without
# XXX - the catastrophic failure.
apt-get install --yes postgresql-client-common postgresql-client
PGPASSWORD=${database_password} psql -c "create database data_refinery" -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}
PGPASSWORD=${database_password} psql -c "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password';" -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}
PGPASSWORD=${database_password} psql -c 'GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user;' -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}
PGPASSWORD=${database_password} psql -c 'ALTER USER data_refinery_user CREATEDB;' -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}
# This won't work on RDS. AWS won't give you the superuser ability.
# PGPASSWORD=${database_password} psql -c 'ALTER ROLE data_refinery_user superuser;' -h ${database_host} -p 5432 -U ${database_user} -d ${database_name}
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


# Create a directory for docker to use as a volume.
mkdir /home/ubuntu/docker_volume
chmod a+rwx /home/ubuntu/docker_volume


# Install Nomad
chmod +x install_nomad.sh
./install_nomad.sh

# Start the Nomad agent in client mode.
nomad agent -config client.hcl > /var/log/nomad_client.log &
