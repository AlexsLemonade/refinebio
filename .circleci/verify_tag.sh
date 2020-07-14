#!/bin/bash

# This script verifies that the tag triggering this deploy was signed
# by a trusted member of the CCDL.

cd ~/refinebio || exit

for key in keys/*; do
    gpg --import "keys/$key"
done

# If it is not a good key then the exit code is 1, which will cause
# the deploy to fail.
git tag --verify "$1"
