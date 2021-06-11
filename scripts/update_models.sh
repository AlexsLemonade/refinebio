#! /bin/sh
# Makes migrations and re-installs so Docker images update locally

# Exit on fail
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of refinebio
cd ..

if ! docker ps | tail -n +2 | awk '{ print $NF }' | grep drdb > /dev/null; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

# Default to "0.0.0.dev" for system version if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="0.0.0.dev$(date +%s)"
    export SYSTEM_VERSION
fi

# Put this in place for common to read from.
echo "$SYSTEM_VERSION" > common/version

# Ensure there is only one distribution to copy over.
rm -f common/dist/*

./scripts/make_migrations.sh && cd common && python3 setup.py sdist
