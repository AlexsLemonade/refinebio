#!/bin/bash

while [ 0 -ne $(PGPASSWORD=data_refinery_password psql \
                          -d data_refinery \
                          -U data_refinery_user \
                          -h 127.0.0.1 \
                          -c 'select 1;' 2>&1 \
                    | grep -c '^ERROR') ]; do
    echo "Postgres isn't alive yet."
    sleep 1
done
