#!/bin/bash

# This script installs Nomad onto an Ubuntu machine. It should be run with `sudo`.
# It also assumes grep, curl, gpg, and sha256sum are installed on the mahcine.
# It will install unzip if it is not already installed.
# These assumptions/installations were based on what is installed by default in
# the AWS Ubuntu AMI.

# Import Hashicorps' Key.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import

# Download the binary and signature files.
curl -0s https://releases.hashicorp.com/nomad/0.9.0-beta2/nomad_0.9.0-beta2_linux_amd64.zip > nomad_0.9.0-beta2_linux_amd64.zip
curl -0s https://releases.hashicorp.com/nomad/0.9.0-beta2/nomad_0.9.0-beta2_SHA256SUMS > nomad_0.9.0-beta2_SHA256SUMS
curl -0s https://releases.hashicorp.com/nomad/0.9.0-beta2/nomad_0.9.0-beta2_SHA256SUMS.sig > nomad_0.9.0-beta2_SHA256SUMS.sig

# Verify the signature file is untampered.
gpg_ok=$(gpg --verify nomad_0.9.0-beta2_SHA256SUMS.sig nomad_0.9.0-beta2_SHA256SUMS |& grep Good)
if [[ "$gpg_ok" = "" ]]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c nomad_0.9.0-beta2_SHA256SUMS |& grep OK)
if [[ "$shasum_ok" = "" ]]; then
    echo "Could not verify the Nomad checksum provided by Hashicorp."
    exit 1
fi

apt-get update
apt-get install -y unzip

unzip -d /usr/bin nomad_0.9.0-beta2_linux_amd64.zip
chmod a+rx /usr/bin/nomad

# Cleanup after ourselves
rm nomad_0.9.0-beta2_linux_amd64.zip
rm nomad_0.9.0-beta2_SHA256SUMS
rm nomad_0.9.0-beta2_SHA256SUMS.sig
