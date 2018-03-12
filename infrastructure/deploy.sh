#!/bin/bash -e

# This script must be run from /infrastructure!

# Make our IP address known to terraform.
source ../common.sh
export TF_VAR_HOST_IP=$(get_ip_address)

# Copy ingress config to top level so it can be applied.
cp deploy/ci_ingress.tf .

# Open up ingress to AWS for Circle, stop jobs, migrate DB.
terraform apply -auto-approve

# Remove the ingress config so the next `terraform apply` will remove
# access for Circle.
rm ci_ingress.tf

# Taint the lead server so that it and all other servers (which depend
# on it) will be torn down and re-initialized.
terraform taint aws_instance.nomad_server_1
terraform apply -auto-approve
