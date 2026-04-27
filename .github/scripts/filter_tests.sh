#!/bin/bash

set -e

git log --format=oneline -n 1 "$GITHUB_SHA"

if [[ $(git log --format=oneline -n 1 "$GITHUB_SHA") = *"noslow"* ]]; then
	echo "Skipping slow tests..."
	./bin/rbio test:workers --exclude-tag=slow "$@"
else
	echo "Running all tests..."
	./bin/rbio test:workers "$@"
fi
