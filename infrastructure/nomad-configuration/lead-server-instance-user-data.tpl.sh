#!/bin/bash

# This is a template for the instance-user-data.sh script for the Lead Nomad Server.
# For more information on instance-user-data.sh scripts, see:
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/user-data.html

# This script will be formatted by Terraform, which will read files from the
# project into terraform variables, and then template them into the following
# script. These will then be written out to files so that they can be used
# locally. This means that any variable referenced as `${name}` is NOT a shell
# variable, it is a template variable for Terraform to fill in. DO NOT treat
# them as normal shell variables.

# A more ideal solution than this would be if we could just give AWS a list of
# files to put onto the instance for us, but they only give us this one script
# to do it with. Nomad has file provisioners which will put files onto the
# instance after it starts up, but those run after this script runs.

# This template varies from nomad-server-instance-user-data.tpl.sh in
# only a few ways.  Because this will be run first, it does not need
# to know the IP of any other server. Instead its IP will be used by
# the other Nomad Servers to join the Raft. Additionally since we only
# need to register the Nomad Jobs once so we do it in this script
# since it will only be run once in total (the other server startup
# script will be run twice, once by each of the other Nomad Servers).


# Change to home directory of the default user
cd /home/ubuntu || exit

# Install, configure and launch our CloudWatch Logs agent
cat <<EOF >awslogs.conf
[general]
state_file = /var/lib/awslogs/agent-state
EOF

mkdir /var/lib/awslogs
wget https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
python ./awslogs-agent-setup.py --region "${region}" --non-interactive --configfile awslogs.conf
# Rotate the logs, delete after 3 days.
echo "
/var/log/nomad_server.log {
    missingok
    notifempty
    compress
    size 20k
    daily
    maxage 3
}" >> /etc/logrotate.conf

# Output the files we need to start up Nomad and register jobs.
# Note that the lines starting with "$" are where
# Terraform will template in the contents of those files.

# Create the script to install Nomad.
cat <<"EOF" > install_nomad.sh
${install_nomad_script}
EOF

# Create the Nomad Server configuration.
cat <<"EOF" > server.hcl
${nomad_server_config}
EOF

# Install Nomad.
chmod +x install_nomad.sh
./install_nomad.sh

# Don't install our graphite/statsd container because it bogs down the
# instance a lot.  After killing graphite on a long-running Nomad
# instance, Nomad HTTP responses went from ~5 seconds to <1 second.
# Based on: https://github.com/graphite-project/docker-graphite-statsd

# Graphite has no auto-deletion, see: https://stackoverflow.com/a/19899267/1135467
# Overwrites defeault set here: https://github.com/graphite-project/docker-graphite-statsd/blob/master/conf/opt/graphite/conf/storage-schemas.conf
# mkdir graphite_configs
# cat <<"EOF" > graphite_configs/storage-schemas.conf
# [everything_max_1d]
# pattern = .*
# retentions = 60s:1d
# EOF

# # Actually run Graphite with our new config
# docker run -d\
#  --name graphite\
#  --restart=always\
#  -v /home/ubuntu/graphite_configs/storage-schemas.conf:/opt/graphite/conf/storage-schemas.conf \
#  -p 80:80\
#  -p 2003-2004:2003-2004\
#  -p 2023-2024:2023-2024\
#  -p 8125:8125/udp\
#  -p 8126:8126\
#  graphiteapp/graphite-statsd

# Start the Nomad agent in server mode via Monit
apt-get -y install monit htop

echo "
#!/bin/bash
nomad status
exit \$?
" >> /home/ubuntu/nomad_status.sh
chmod +x /home/ubuntu/nomad_status.sh

echo "
#!/bin/bash
killall nomad
sleep 120
nomad agent -config /home/ubuntu/server.hcl > /var/log/nomad_server.log &
" >> /home/ubuntu/kill_restart_nomad.sh
chmod +x /home/ubuntu/kill_restart_nomad.sh
/home/ubuntu/kill_restart_nomad.sh

echo '
check program nomad with path "/bin/bash /home/ubuntu/nomad_status.sh" as uid 0 and with gid 0
    start program = "/bin/bash /home/ubuntu/kill_restart_nomad.sh" as uid 0 and with gid 0 with timeout 240 seconds
    if status != 0
        then restart
set daemon 300
' >> /etc/monit/monitrc

service monit restart

# Create the CW metric job in a crontab
# write out current crontab
crontab -l > tempcron
# shellcheck disable=2016
echo -e '#!/bin/bash\naws cloudwatch put-metric-data --metric-name NomadQueueLength --namespace ${user}-${stage} --value `nomad status | grep dispatch | grep -e pending -e running | wc -l` --region ${region}' >> update_metric.sh
chmod +x update_metric.sh

# echo new cron into cron file
echo -e "SHELL=/bin/bash\nPATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin\n* * * * * /home/ubuntu/update_metric.sh" >> tempcron
# install new cron file
crontab tempcron
rm tempcron
