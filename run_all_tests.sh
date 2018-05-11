#!/bin/bash

# Runs all common and app specific tests.
# End to end tests require Nomad to be running for the test environment.
# This can be done with:
# sudo -E ./run_nomad.sh -e test
# This should mirror what happens in the CircleCI config (.circleci/config).
./api/run_tests.sh
./common/run_tests.sh
./foreman/run_tests.sh
./workers/run_tests.sh
