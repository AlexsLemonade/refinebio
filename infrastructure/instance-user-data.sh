#!/bin/bash
set -e

echo 'DOCKER_STORAGE_OPTIONS="--storage-driver overlay2"' > /etc/sysconfig/docker-storage
/etc/init.d/docker restart

# This is the public key from above - one-time step.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import

# Download the binary and signature files.
wget https://releases.hashicorp.com/nomad/0.7.1/nomad_0.7.1_linux_amd64.zip
wget https://releases.hashicorp.com/nomad/0.7.1/nomad_0.7.1_SHA256SUMS
wget https://releases.hashicorp.com/nomad/0.7.1/nomad_0.7.1_SHA256SUMS.sig

# Verify the signature file is untampered.
if [ $(gpg --verify nomad_0.7.1_SHA256SUMS.sig nomad_0.7.1_SHA256SUMS 2>&1 grep Good) = "" ]; then
    >&2 echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
if [ $(shasum -a 256 -c nomad_0.7.1_SHA256SUMS 2>&1 | grep OK) = ""]; then
    >&2 echo "Could not verify the Nomad checksum provided by Hashicorp."
    exit 1
fi

unzip -d /usr/bin nomad_0.7.1_linux_amd64.zip
