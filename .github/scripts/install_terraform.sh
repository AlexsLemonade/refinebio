#!/bin/bash

TERRAFORM_VERSION=0.11.8
curl -0s "https://releases.hashicorp.com/terraform/${TERRAFORM_VERSION}/terraform_${TERRAFORM_VERSION}_linux_amd64.zip" \
    > "terraform_${TERRAFORM_VERSION}_linux_amd64.zip"
sudo unzip -d /usr/bin "terraform_${TERRAFORM_VERSION}_linux_amd64.zip"
sudo chmod a+rx /usr/bin/terraform
