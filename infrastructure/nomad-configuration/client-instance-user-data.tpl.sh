MIME-Version: 1.0
Content-Type: multipart/mixed; boundary="==MYBOUNDARY=="

--==MYBOUNDARY==
Content-Type: text/x-shellscript; charset="us-ascii"

# This is a template for the instance-user-data.sh script for the nomad client.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from the
# project into terraform variables, and then template them into the following
# script. These will then be written out to files so that they can be used
# locally. This means that any variable referenced using `{}` is NOT a shell
# variable, it is a template variable for Terraform to fill in. DO NOT treat
# them as normal shell variables.

# A more ideal solution than this would be if we could just give AWS a list of
# files to put onto the instance for us, but they only give us this one script
# to do it with. Nomad has file provisioners which will put files onto the
# instance after it starts up, but those run after this script runs.

yum -y update
yum install -y jq iotop dstat speedometer awscli docker.io htop wget

ulimit -n 65536

# Configure and mount the EBS volume
mkdir -p /var/ebs/

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
