### AMI Creation

This directory contains the terraform config for an instance "AMI Template
Instance" which upgrades Ubuntu and installs docker and nomad.
This can be run with `terraform apply`.

You can then use the script `create_ubuntu_ami.sh` to create a new Ubuntu AMI with the name
`ccdl-ubuntu-18.04-<TIMESTAMP>`.
This is used by the foreman, the API, and pg_bouncer.

Or you can use the script `create_ecs_ami.sh` to create a new [ECS-optimized](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-optimized_AMI.html) AMI with the name
`ccdl-ubuntu-18.04-<TIMESTAMP>`.
This is used by worker instances of the AWS Batch cluster.

Be sure that everything on the template instance is fully installed before using either of the
create AMI scripts, and that you only destroy terraform until the AMI is complete.
