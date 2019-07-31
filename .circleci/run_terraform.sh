#!/bin/bash -e

# Import Hashicorps' Key.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import


# Install terraform and nomad
cd
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

sudo apt-get update
sudo apt-get install lxc -y  # Install lxc, which is required by nomad

NOMAD_VERSION=v0.9.3
wget -N https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_linux_amd64.zip
wget -N https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_SHA256SUMS
wget -N https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_SHA256SUMS.sig


# Verify the signature file is untampered.
gpg_ok=$(gpg --verify nomad_${NOMAD_VERSION}_SHA256SUMS.sig nomad_${NOMAD_VERSION}_SHA256SUMS |& grep Good)
if [[ "$gpg_ok" == "" ]]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c nomad_${NOMAD_VERSION}_SHA256SUMS |& grep OK)
if [[ "$shasum_ok" == "" ]]; then
    echo "Could not verify the Nomad checksum provided by Hashicorp."
    exit 1
fi

unzip nomad_${NOMAD_VERSION}_linux_amd64.zip
sudo mv nomad /usr/local/bin/

cd ~/refinebio/infrastructure

# Circle won't set the branch name for us, so do it ourselves.
source ~/refinebio/common.sh
branch=$(get_master_or_dev $CIRCLE_TAG)

if [[ "$CIRCLE_TAG" == *"-hotfix" && $branch == "dev" ]]; then
    ENVIRONMENT=prod
elif [[ $branch == "master" ]]; then
    ENVIRONMENT=prod
elif [[ $branch == "dev" ]]; then
    ENVIRONMENT=staging
else
    echo "Why in the world was run_terraform.sh called from a branch other than dev or master?!?!?"
    exit 1
fi

# New deployment
./deploy.sh -e $ENVIRONMENT -v $CIRCLE_TAG -u circleci
