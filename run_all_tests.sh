#!/bin/bash

# Runs all common and app specific tests.
# End to end tests require Nomad to be running for the test environment.
# This can be done with:
# sudo -E ./run_nomad.sh -e test
# This should mirror what happens in the CircleCI config (.circleci/config).

# Ensure that Nomad is running first
if ps aux | grep test_nomad | grep -v grep; then
    mkdir -p test_volume && chmod -R a+rw test_volume
    ./update_models.sh
    ./api/run_tests.sh
    ./common/run_tests.sh
    ./foreman/run_tests.sh
    ./workers/run_tests.sh
else
    echo "You must start the nomad test environment first with" >&2
    echo "'sudo -E ./run_nomad.sh -e test'" >&2
    exit 1
fi