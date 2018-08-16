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

apt-get -y update
sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt/ $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key -y add -
apt-get -y update
apt-get -y install pgbouncer postgresql-client --allow-unauthenticated

# Set up PG Bouncer
cat <<"EOF" > /etc/pgbouncer/pgbouncer.ini
[databases]
${database_name} = host=${database_host} port=5430 dbname=${database_name}

[pgbouncer]
listen_addr = *
listen_port = 5432
auth_type = md5
auth_file = /etc/pgbouncer/userlist.txt
pool_mode = transaction
server_reset_query =
pidfile=/tmp/pgbouncer.pidfile
EOF

# Set up PG Bouncer
cat <<"EOF" > /etc/pgbouncer/userlist.txt
"${database_user}" "${database_password}"
EOF

echo "ulimit -n 16384" >> /etc/default/pgbouncer

sudo service pgbouncer start

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi