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
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get upgrade -y

# Set up the docker apt repository
apt-get install -y apt-transport-https ca-certificates curl gnupg-agent software-properties-common

apt-key add - <<EOF
${docker_apt_key}
EOF

# Adding the PPA for python3.5
add-apt-repository ppa:deadsnakes/ppa

add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
apt-get update -y

# This is the latest version of python that the aws logs can use (don't work with 3.7)
apt-get install --yes python3.5

# Install our dependencies
apt-get install --yes jq iotop dstat speedometer awscli docker-ce docker-ce-cli containerd.io chrony htop monit

# Add the regular user to the docker group so they can perform docker commands
sudo usermod -aG docker ubuntu

# Change to home directory of the default user
cd /home/ubuntu || exit

# Set up the AWS NTP
# via https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/set-time.html#configure_ntp
echo 'server 169.254.169.123 prefer iburst' | cat - /etc/chrony/chrony.conf > temp && mv temp /etc/chrony/chrony.conf
/etc/init.d/chrony restart

# We seem to have run into: https://medium.com/spaceapetech/what-the-arp-is-going-on-b4bc0e73e4d4
echo 'net.ipv4.neigh.default.gc_thresh1 = 0' | tee /etc/sysctl.d/55-arp-gc_thresh1.conf

# Restart to apply the updates
reboot now
