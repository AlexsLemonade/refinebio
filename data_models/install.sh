#!/bin/bash

# Run this script with 'sudo -u postgres install.sh'

# This script defaults to postgres version 9.5.
# To use a different version set the POSTGRES_VERSION env var
# (i.e. 'sudo POSTGRES_VERSION=9.4 -u postgres install.sh')
if [ -z $POSTGRES_VERSION ]; then POSTGRES_VERSION=9.5; fi

createdb bioinformatics_mill

psql bioinformatics_mill -c "CREATE ROLE bioinformatics_mill_user WITH LOGIN PASSWORD 'bioinformatics_password';"
psql bioinformatics_mill -c 'GRANT ALL PRIVILEGES ON DATABASE bioinformatics_mill TO bioinformatics_mill_user;'
psql bioinformatics_mill -c 'ALTER USER bioinformatics_mill_user CREATEDB;'

iptables -A INPUT -s 172.17.0.0/16 -j ACCEPT

echo 'host    all             all             172.17.0.0/16           md5' >> /etc/postgresql/$POST_GRES_VERSION/main/pg_hba.conf

sed -i "s/#\(listen_addresses\)\( = 'localhost'\)/\1 = '*'/" /etc/postgresql/$POSTGRES_VESRION/main/postgresql.conf
