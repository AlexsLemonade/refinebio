#!/bin/bash

# Exit on failure.
set -e

git log --format=oneline -n 1 "$GITHUB_SHA"
if [[ $(git log --format=oneline -n 1 "$GITHUB_SHA") = *"noslow"* ]]; then
	echo "Skipping slow tests..."
	./workers/run_tests.sh --exclude-tag=slow "$@"
else
	echo "Running all tests..."
	./workers/run_tests.sh "$@"
fi

./.github/scripts/upload_test_coverage.sh workers
