# This script can be used as the instance-user-data.sh script for AWS
# instances to configure a new instance to be used as a deploy box.

# To create a new deploy box (there only ever needs to be one at a
# time.), there are few steps to follow:
#  - Make sure there is a deploy VPC. If there isn't create one so the deploy box is isolated.
#  - Make sure the deploy VPC has a security group. If it doesn't create one and allow ingress
#    TCP connections on port 22 from 0.0.0.0/0
#  - Create a new EC2 instance with the following settings:
#    - Select an Ubuntu 18.04 HVM SSD AMI for the region. If you are using us-east-1 then
#      ubuntu/images/hvm-ssd/ubuntu-bionic-18.04-amd64-server-20180912 (ami-0ac019f4fcb7cb7e6)
#      should be correct.
#    - A t3.small instance should be sufficient.
#    - Set the EBS volume to have at least 100 GB.
#    - Use the deploy VPC and security group referenced in earlier steps.
#    - Paste the contents of this script into the instance-user-data.

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

apt-key fingerprint 9DC858229FC7DD38854AE2D88D81803C0EBFCD88

add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

# Install git-crypt
apt-get update -qq
apt-get install -y \
        make \
        g++ \
        libssl-dev \
        python3 \
        python3-pip \
        jq \
        docker-ce \
        unzip

groupadd docker
usermod -aG docker ubuntu

cd /home/ubuntu
git clone https://github.com/AGWA/git-crypt.git
cd git-crypt
make
make install

touch /var/log/docker_update.log
chown ubuntu:ubuntu /var/log/docker_update.log

touch /var/log/deploy.log
chown ubuntu:ubuntu /var/log/deploy.log

# Checkout the repo onto the box.
cd /home/ubuntu
git clone https://github.com/AlexsLemonade/refinebio.git
