#!/bin/bash

# Runs all common and app specific tests.
# This should mirror what happens in the CircleCI config (.circleci/config).
./common/run_tests.sh
./foreman/run_tests.sh
./workers/run_tests.sh

