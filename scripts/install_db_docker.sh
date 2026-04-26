#! /bin/sh

# Bootstrap the data_refinery database in the running drdb container.
# Postgres trusts local socket connections, so no PGPASSWORD needed.

set -e

docker exec drdb psql -U postgres -c "CREATE DATABASE data_refinery"
docker exec drdb psql -U postgres -c "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password'"
docker exec drdb psql -U postgres -c 'GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user'
docker exec drdb psql -U postgres -c 'ALTER USER data_refinery_user CREATEDB'
docker exec drdb psql -U postgres -c 'ALTER ROLE data_refinery_user superuser'
docker exec drdb psql -U postgres -d data_refinery -c 'CREATE EXTENSION IF NOT EXISTS hstore'
