#! /bin/bash

docker run -it -e PGPASSWORD=mysecretpassword --rm --link drdb:postgres postgres:9.6.6 psql -c "create database data_refinery" -h postgres -U postgres
docker run -it -e PGPASSWORD=mysecretpassword --rm --link drdb:postgres postgres:9.6.6 psql -c "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password';" -h postgres -U postgres
docker run -it -e PGPASSWORD=mysecretpassword --rm --link drdb:postgres postgres:9.6.6 psql -c 'GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user;' -h postgres -U postgres
docker run -it -e PGPASSWORD=mysecretpassword --rm --link drdb:postgres postgres:9.6.6 psql -c 'ALTER USER data_refinery_user CREATEDB;' -h postgres -U postgres
