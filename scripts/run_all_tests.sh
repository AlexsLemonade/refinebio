#!/bin/sh

# Runs all common and app specific tests.
# End to end tests require Nomad to be running for the test environment.
# This can be done with:
# sudo -E ./scripts/run_nomad.sh -e test
# This should mirror what happens in the CircleCI config (.circleci/config).

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of the refinebio project
cd ..

# Ensure that Nomad is running first
# The -f option in `pgrep` means that pgrep checks the full command line
if ! pgrep -f test_nomad > /dev/null; then
    echo "You must start the nomad test environment first with" >&2
    echo "sudo -E ./scripts/run_nomad.sh -e test" >&2
    exit 1
# Then ensure postgres is running
elif ! docker ps | tail -n +2 | awk '{ print $NF }' | grep drdb > /dev/null; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

mkdir -p test_volume && chmod -R a+rw test_volume
./scripts/update_models.sh
./api/run_tests.sh
./common/run_tests.sh
./foreman/run_tests.sh
./workers/run_tests.sh
