#!/bin/bash -e

if [ -z $TF_VAR_user ]; then
    echo 'You must set the $TF_VAR_user environment variable!' >&2
    exit 1
fi

if [ -z $TF_VAR_stage ]; then
    STAGE=dev
else
    STAGE=$TF_VAR_stage
fi

# Terraform leaves some stack-specific config lying around.
# This config messes up initialization of an arbitrary stage, since it is stack-specific.
# Therefore in order to be able to switch between different stacks, just remove that config.
rm -rf .terraform

# Force copy will allow us to get past the prompt asking if we want to
# switch to the remote backend.
terraform init \
    -force-copy \
    -backend-config="bucket=refinebio-tfstate-deploy-$STAGE" \
    -backend-config="key=terraform-${TF_VAR_user}.tfstate" \
    -backend-config="dynamodb_table=refinebio-terraform-lock"
