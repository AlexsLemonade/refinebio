#!/bin/bash

# Exit on failure
set -e

if [[ $(git log --format=oneline -n 1 "$GITHUB_SHA") = *"noslow"* ]]; then
	echo "Skipping slow tests..";
	./workers/run_tests.sh --exclude-tag=slow "$@"
else
	echo "Running all tests..";
	./workers/run_tests.sh "$@"
fi

./.circleci/upload_test_coverage.sh workers
