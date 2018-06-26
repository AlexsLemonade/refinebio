#! /bin/bash -e
# Makes migrations and re-installs so Docker images update locally

if ! docker ps | tail -n +2 | awk '{ print $NF }' | grep drdb > /dev/null; then
    echo "You must start Postgres first with './run_postgres.sh'"
    exit 1
fi

./common/make_migrations.sh && cd common && python setup.py sdist
