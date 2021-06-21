#!/bin/sh

# This script installs Nomad onto an Ubuntu machine. It should be run with `sudo`.
# It also assumes grep, curl, gpg, and sha256sum are installed on the mahcine.
# It will install unzip if it is not already installed.
# These assumptions/installations were based on what is installed by default in
# the AWS Ubuntu AMI.

# First, update the local CA certificates
# adapted from: https://stackoverflow.com/questions/62107565/wget-and-curl-stopped-working-with-https-wrongly-complain-about-an-expired-cert
apt update && apt install ca-certificates
sed -i '/mozilla\/AddTrust_External/d' /etc/ca-certificates.conf
update-ca-certificates -f -v

# Import Hashicorps' Key.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import

# Download the binary and signature files.
RELEASE_VERSION=0.11.3
RELEASE_PLATFORM="$(uname | tr '[:upper:]' '[:lower:]')"
export RELEASE_VERSION RELEASE_PLATFORM

sha256sum() {
    if [ "$RELEASE_PLATFORM" = "darwin" ]; then
        shasum -a 256 "$@" ;
    else
        shasum "$@"
    fi
}

curl -0s "https://releases.hashicorp.com/nomad/${RELEASE_VERSION}/nomad_${RELEASE_VERSION}_${RELEASE_PLATFORM}_amd64.zip" > "nomad_${RELEASE_VERSION}_${RELEASE_PLATFORM}_amd64.zip"
curl -0s "https://releases.hashicorp.com/nomad/${RELEASE_VERSION}/nomad_${RELEASE_VERSION}_SHA256SUMS" > "nomad_${RELEASE_VERSION}_SHA256SUMS"
curl -0s "https://releases.hashicorp.com/nomad/${RELEASE_VERSION}/nomad_${RELEASE_VERSION}_SHA256SUMS.sig" > "nomad_${RELEASE_VERSION}_SHA256SUMS.sig"

# Verify the signature file is untampered.
gpg_ok=$(gpg --verify "nomad_${RELEASE_VERSION}_SHA256SUMS.sig" "nomad_${RELEASE_VERSION}_SHA256SUMS" 2>&1 | grep Good)
if [ "$gpg_ok" = "" ]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c "nomad_${RELEASE_VERSION}_SHA256SUMS" 2>&1 | grep OK)
if [ "$shasum_ok" = "" ]; then
    echo "Could not verify the Nomad checksum provided by Hashicorp."
    exit 1
fi

if [ "$RELEASE_PLATFORM" = "linux" ]; then
    apt-get update
    apt-get install -y unzip
    unzip -d /usr/bin "nomad_${RELEASE_VERSION}_${RELEASE_PLATFORM}_amd64.zip"
    chmod a+rx /usr/bin/nomad
fi
if [ "$RELEASE_PLATFORM" = "darwin" ]; then
    unzip -d /usr/local/bin "nomad_${RELEASE_VERSION}_${RELEASE_PLATFORM}_amd64.zip"
    chmod a+rx /usr/local/bin/nomad
fi

# Cleanup after ourselves
rm "nomad_${RELEASE_VERSION}_${RELEASE_PLATFORM}_amd64.zip"
rm "nomad_${RELEASE_VERSION}_SHA256SUMS"
rm "nomad_${RELEASE_VERSION}_SHA256SUMS.sig"
