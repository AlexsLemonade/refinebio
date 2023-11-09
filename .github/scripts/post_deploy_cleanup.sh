#!/bin/bash

set -e

# TODO: this should also be responsible for running another
# deploy to remove ci_ingress after the end to end tests run.

# Don't leave secrets lying around.
# shellcheck disable=SC2029
ssh -o StrictHostKeyChecking=no \
    -o ServerAliveInterval=15 \
    -i infrastructure/data-refinery-key.pem \
    "ubuntu@${DEPLOY_BOX_IP}" \
    "cd refinebio && git clean -f"
