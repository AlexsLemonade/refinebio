#! /bin/sh

docker run \
    --env PGPASSWORD=mysecretpassword \
    --interactive \
    --link drdb:postgres \
    --rm \
    --tty \
    postgres:9.6.6 \
    psql -h postgres -U postgres -d data_refinery
