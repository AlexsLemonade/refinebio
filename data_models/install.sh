#!/bin/bash

# Run this script with 'sudo -u postgres install.sh'

# This script defaults to postgres version 9.5.
# To use a different version set the POSTGRES_VERSION env var
# (i.e. 'sudo POSTGRES_VERSION=9.4 -u postgres install.sh')
if [ -z $POSTGRES_VERSION ]; then POSTGRES_VERSION=9.5; fi

createdb data_refinery

psql data_refinery -c "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password';"
psql data_refinery -c 'GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user;'
psql data_refinery -c 'ALTER USER data_refinery_user CREATEDB;'

# See https://www.howtogeek.com/177621/the-beginners-guide-to-iptables-the-linux-firewall/
# for more information.
iptables -A INPUT -s 172.17.0.0/16 -j ACCEPT

# See http://stackoverflow.com/questions/31249112/allow-docker-container-to-connect-to-a-local-host-postgres-database
# and https://blog.jsinh.in/how-to-enable-remote-access-to-postgresql-database-server/#.WNFiBCErLmG
# for more information on these two settings.
echo 'host    all             all             172.17.0.0/16           md5' >> /etc/postgresql/$POST_GRES_VERSION/main/pg_hba.conf

sed -i "s/#\(listen_addresses\)\( = 'localhost'\)/\1 = '*'/" /etc/postgresql/$POSTGRES_VESRION/main/postgresql.conf
