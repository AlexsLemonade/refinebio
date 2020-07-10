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

add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
apt-get update -y

# Install our dependencies
apt-get install --yes jq iotop dstat speedometer awscli docker-ce docker-ce-cli containerd.io chrony htop monit

# Add the regular user to the docker group so they can perform docker commands
sudo usermod -aG docker ubuntu

# Change to home directory of the default user
cd /home/ubuntu || exit

# Install nomad
cat <<"EOF" > install_nomad.sh
${install_nomad_script}
EOF
chmod +x install_nomad.sh
./install_nomad.sh

# Set up the AWS NTP
# via https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/set-time.html#configure_ntp
echo 'server 169.254.169.123 prefer iburst' | cat - /etc/chrony/chrony.conf > temp && mv temp /etc/chrony/chrony.conf
/etc/init.d/chrony restart

# Restart to apply the updates
reboot now
