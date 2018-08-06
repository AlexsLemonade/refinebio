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

# Force copy will allow us to get past the prompt asking if we want to
# switch to the remote backend.
terraform init \
    -force-copy \
    -backend-config="bucket=refinebio-tfstate-deploy-$STAGE" \
    -backend-config="key=terraform-${TF_VAR_user}.tfstate"
