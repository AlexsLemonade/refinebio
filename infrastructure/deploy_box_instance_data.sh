#!/bin/bash -e

# This script can be used as the instance-user-data.sh script for AWS
# instances to configure a new instance to be used as a deploy box.


# To create a new deploy box (there only ever needs to be one at a
# time.), there are few steps to follow using the AWS console:
#  - Make sure there is a deploy VPC. If there isn't create one so the deploy box is isolated.
#    See https://docs.aws.amazon.com/directoryservice/latest/admin-guide/gsg_create_vpc.html for detailed instructions.

#  - Make sure the deploy VPC has a security group. If it doesn't create one and allow ingress
#    TCP connections on port 22 from 0.0.0.0/0
#    See https://docs.aws.amazon.com/vpc/latest/userguide/VPC_SecurityGroups.html for detailed instructions.

#  - Make an AWS Key Pair named `data-refinery-deployer` by importing the data-refinery-key.pem
#    key (which is in the same directory as this script)into AWS if it does not already exist:
#    See https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#how-to-generate-your-own-key-and-import-it-to-aws
#    for detailed instructions.

#  - Create a new EC2 instance with the following settings:
#    - Select an Ubuntu 18.04 HVM SSD AMI for the region. If you are using us-east-1 then
#      ubuntu/images/hvm-ssd/ubuntu-bionic-18.04-amd64-server-20180912 (ami-0ac019f4fcb7cb7e6)
#      should be correct.
#    - A t3.medium instance should be sufficient.
#    - Set the EBS volume to have at least 100 GB.
#    - Use the deploy VPC and security group referenced in earlier steps.
#    - It needs to have a public IP address (so set 'Auto-assign Public IP' to 'Enable').
#    - Paste the contents of this script into the instance-user-data.
#    - When it prompts you for a key, specify the `data-refinery-deployer` key you created earlier.
#    See https://docs.aws.amazon.com/efs/latest/ug/gs-step-one-create-ec2-resources.html for detailed instructions.

#  - Record the IPv4 Public IP of the instance you just created.


# Finally, go into the CirlceCI web application and select the refinebio project.
# Go to the project settings and navigate to the Environment Variables tab.
# Click 'Add Variable' and set the name to DEPLOY_IP_ADDRESS and the value to the
# IP address of the EC2 instance you created.

# Also, if you want to be notified on slack after the deploy finishes, you can add
# another environment variable named ENGAGEMENTBOT_WEBHOOK with the url to a slack webhook.

# Now any deploys that are triggered will run on the instance you created.

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
        unzip \
        postgresql-client

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
chown -R ubuntu:ubuntu refinebio
