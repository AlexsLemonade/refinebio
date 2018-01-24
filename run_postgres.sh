#! /bin/bash

# via https://hub.docker.com/_/postgres/
# 9.6.6 is the current (as of Jan 23 2018) RDS most recent version.
# Password can be exposed to git/CI because this is only for dev/testing purposes, not real data.
docker run --name drdb -e POSTGRES_PASSWORD=mysecretpassword -d postgres:9.6.6
