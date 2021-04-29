#!/bin/bash -e

# Don't leave secrets lying around.
# shellcheck disable=SC2029
ssh -o StrictHostKeyChecking=no \
    -o ServerAliveInterval=15 \
    -i infrastructure/data-refinery-key.pem \
    ubuntu@"${DEPLOY_IP_ADDRESS}" "cd refinebio && git clean -f"
