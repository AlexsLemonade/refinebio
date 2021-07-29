MIME-Version: 1.0
Content-Type: multipart/mixed; boundary="==MYBOUNDARY=="

--==MYBOUNDARY==
Content-Type: text/x-shellscript; charset="us-ascii"

# This is a template for the instance-user-data.sh script for the worker instances.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform.

yum -y update
yum install -y jq iotop dstat speedometer awscli docker.io htop wget

ulimit -n 65536

# Add docker login to the ECS agent configuration. This is based off of
# https://aws.amazon.com/blogs/aws/ec2-container-service-ecs-update-access-private-docker-repos-mount-volumes-in-containers/,
# although it also changes ECS_CONTAINER_STOP_TIMEOUT based on
# https://github.com/aws/amazon-ecs-agent/issues/2312.
# systemctl stop ecs

aws s3 cp s3://data-refinery-secrets/instance-ecs-agent.config instance-ecs-agent.config
cat instance-ecs-agent.config >> /etc/ecs/ecs.config
rm instance-ecs-agent.config

systemctl restart docker --no-block

# Configure and mount the EBS volume
mkdir -p /var/ebs/

# We want to mount the biggest volume that its attached to the instance
# The size of this volume can be controlled with the varialbe
# `volume_size_in_gb` from the file `variables.tf`
ATTACHED_AS=`lsblk -n --sort SIZE | tail -1 | cut -d' ' -f1`

# grep -v ext4: make sure the disk is not already formatted.
if file -s /dev/$ATTACHED_AS | grep data | grep -v ext4; then
	mkfs -t ext4 /dev/$ATTACHED_AS # This is slow
fi
mount /dev/$ATTACHED_AS /var/ebs/

export STAGE="${stage}"
export USER="${user}"
EBS_VOLUME_INDEX="$(wget -q -O - http://169.254.169.254/latest/meta-data/instance-id)"

chown ec2-user:ec2-user /var/ebs/
echo "$EBS_VOLUME_INDEX" > /var/ebs/VOLUME_INDEX
chown ec2-user:ec2-user /var/ebs/VOLUME_INDEX

# Set up the required database extensions.
# HStore allows us to treat object annotations as pseudo-NoSQL data tables.
amazon-linux-extras install postgresql10
PGPASSWORD=${database_password} psql -c 'CREATE EXTENSION IF NOT EXISTS hstore;' -h "${database_host}" -p 5432 -U "${database_user}" -d "${database_name}"

# Change to home directory of the default user
cd /home/ec2-user || exit

# Delete the cloudinit and syslog in production.
if [[ $STAGE = *"prod"* ]]; then
    rm /var/log/cloud-init.log
    rm /var/log/cloud-init-output.log
    rm /var/log/syslog
fi




--==MYBOUNDARY==--
