### AMI Creation

This directory contains the terraform config for an instance "AMI Template
Instance" which upgrades Ubuntu and installs docker and nomad. You can then use
the script `create_ami.sh` to create a new AMI with the name
`ccdl-ubuntu-18.04-<TIMESTAMP>`.
