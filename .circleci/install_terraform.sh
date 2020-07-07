#!/bin/bash

cd || exit
TERRAFORM_VERSION=0.11.8
wget -N https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_linux_amd64.zip
wget -N https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_SHA256SUMS
wget -N https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_SHA256SUMS.sig


# Verify the signature file is untampered.
gpg_ok=$(gpg --verify terraform_${TERRAFORM_VERSION}_SHA256SUMS.sig terraform_${TERRAFORM_VERSION}_SHA256SUMS |& grep Good)
if [[ "$gpg_ok" == "" ]]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c terraform_${TERRAFORM_VERSION}_SHA256SUMS |& grep OK)
if [[ "$shasum_ok" == "" ]]; then
    echo "Could not verify the Terraform checksum provided by Hashicorp."
    exit 1
fi

unzip terraform_${TERRAFORM_VERSION}_linux_amd64.zip
sudo mv terraform /usr/local/bin/
