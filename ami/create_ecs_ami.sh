#!/bin/sh

instance_id="$(aws ec2 describe-instances | jq -r -f instance_id.jq)"

if [ -z "$instance_id" ]; then
    echo "Could not find an instance with name 'AMI Template Instance'"
fi

echo "Creating an ami for instance $instance_id..."
ami_name="ccdl-ecs-optimized-$(date "+%Y-%m-%dT%H.%M.%S")"
aws ec2 create-image --instance-id "$instance_id" --name "$ami_name"
echo "Created an ami with name $ami_name"
