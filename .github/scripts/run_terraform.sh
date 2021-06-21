#!/bin/bash -e

# Import Hashicorps' Key.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import


# Install terraform and nomad
cd
TERRAFORM_VERSION=0.13.5
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

cd ~/refinebio/infrastructure

# Circle won't set the branch name for us, so do it ourselves.
source ~/refinebio/scripts/common.sh
branch=$(get_master_or_dev "$CI_TAG")

if [[ $branch == "master" ]]; then
    ENVIRONMENT=prod
    BATCH_USE_ON_DEMAND_INSTANCES="false"
elif [[ $branch == "dev" ]]; then
    ENVIRONMENT=staging
    BATCH_USE_ON_DEMAND_INSTANCES="true"
else
    echo "Why in the world was run_terraform.sh called from a branch other than dev or master?!?!?"
    exit 1
fi

# New deployment (use -u circleci since we used to run on CircleCI and we don't
# want to recreate all of our resources)
./deploy.sh -e "$ENVIRONMENT" -v "$CI_TAG" -u circleci -i "$BATCH_USE_ON_DEMAND_INSTANCES"
