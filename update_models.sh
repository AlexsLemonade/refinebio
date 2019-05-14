#! /bin/bash -e
# Makes migrations and re-installs so Docker images update locally

if ! docker ps | tail -n +2 | awk '{ print $NF }' | grep drdb > /dev/null; then
    echo "You must start Postgres first with:" >&2
    echo "./run_postgres.sh" >&2
    exit 1
fi

# Default to "local" for system version if we're not running in the cloud.
if [[ -z "$SYSTEM_VERSION" ]]; then
    export SYSTEM_VERSION="local$(date +%s)"
fi

# Put this in place for common to read from.
echo "$SYSTEM_VERSION" > common/version

# Ensure there is only one distribution to copy over.
rm -f common/dist/*

./common/make_migrations.sh && cd common && python setup.py sdist
