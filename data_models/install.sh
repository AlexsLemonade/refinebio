#!/bin/bash

# Run this script with 'sudo -u postgres install.sh'

createdb bioinformatics_mill

psql bioinformatics_mill -c "CREATE ROLE bioinformatics_mill_user WITH LOGIN PASSWORD 'bioinformatics_password';"
psql bioinformatics_mill -c 'GRANT ALL PRIVILEGES ON DATABASE bioinformatics_mill TO bioinformatics_mill_user;'
psql bioinformatics_mill -c 'ALTER USER bioinformatics_mill_user CREATEDB;'
