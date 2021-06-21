#!/bin/bash

# This is a template for the instance-user-data.sh script for the Nomad Server.
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

# Change to home directory of the default user
cd /home/ubuntu || exit

service postgresql stop

apt-get -y update
apt-get -y install pgbouncer postgresql-client

# Set up PG Bouncer
# The values of max_client_conn and default_pool_size.
# This should be a function of the RDS instance max_connections,
# which is determined by instance type.
# Ex., m4.xlarge: 1320 - (.2 * 1320) = 1056
# So, for a pool of 10, 1056 * 10 = 10560
# However, this is basically a dark art and we'll probably want to tweak in the future.

cat << FOE >> /etc/pgbouncer/pgbouncer.ini
[databases]
* = host=${database_host} port=${database_port} user=${database_user} password=${database_password} dbname=${database_name} client_encoding=UNICODE

[pgbouncer]
listen_addr = *
max_client_conn = 10560
default_pool_size = 10
listen_port = ${listen_port}
auth_type = trust
auth_file = /etc/pgbouncer/userlist.txt
pool_mode = transaction
server_reset_query = DISCARD ALL
syslog = 1
pidfile = /var/run/postgresql/pgbouncer.pid
unix_socket_dir = /var/run/postgresql
FOE

# Set up PG Bouncer
# TODO: Actually make the password required to connect through PGBouncer.
cat << FOE >> /etc/pgbouncer/userlist.txt
"${database_user}" "${database_password}"
FOE

echo "ulimit -n 16384" >> /etc/default/pgbouncer


service pgbouncer stop
service pgbouncer start

# Delete the cloudinit and syslog in production.
export STAGE=${stage}
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi
