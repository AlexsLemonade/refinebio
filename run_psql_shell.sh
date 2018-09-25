#! /bin/bash
docker run -it -e PGPASSWORD=mysecretpassword --rm --link drdb:postgres postgres:9.6.6 psql -h postgres -U postgres
