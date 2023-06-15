#! /bin/sh

POSTGRES_VERSION="9.6.6"

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c "create database data_refinery" -h postgres -U postgres

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password';" -h postgres -U postgres

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c 'GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user;' -h postgres -U postgres

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c 'ALTER USER data_refinery_user CREATEDB;' -h postgres -U postgres

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c 'ALTER ROLE data_refinery_user superuser;' -h postgres -U postgres

docker run \
    --env PGPASSWORD=mysecretpassword \
    --link drdb:postgres \
    --rm \
    "postgres:$POSTGRES_VERSION" \
    psql -c 'CREATE EXTENSION IF NOT EXISTS hstore;' -h postgres -U postgres -d data_refinery
