#!/bin/bash -x

# This is a template for the instance-user-data.sh script for the AMI template.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from the
# project into terraform variables, and then template them into the following
# script. These will then be written out to files so that they can be used
# locally. This means that any variable referenced using `{}` is NOT a shell
# variable, it is a template variable for Terraform to fill in. DO NOT treat
# them as normal shell variables.

# First update the base system
yum upgrade -y

# Install our dependencies
sudo amazon-linux-extras install python3

## Required for monit
yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm

yum install -y awscli dstat htop iotop jq monit

# Add the regular user to the docker group so they can perform docker commands
sudo usermod -aG docker ec2-user

# Change to home directory of the default user
cd /home/ec2-user || exit

# We seem to have run into: https://medium.com/spaceapetech/what-the-arp-is-going-on-b4bc0e73e4d4
echo 'net.ipv4.neigh.default.gc_thresh1 = 0' | tee /etc/sysctl.d/55-arp-gc_thresh1.conf

# Restart to apply the updates
reboot now
