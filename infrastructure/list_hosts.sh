#! /bin/bash
# Use like:
# ./list_hosts rjones > hosts
# ./connect_cluster.sh

aws ec2 describe-instances --filters "Name=tag:User,Values=$1" | jq '.Reservations[0].Instances' | grep PublicDnsName | tr -d '"' | sed "s/PublicDnsName: //g" | tr -d "," | awk '{$1=$1};1' | uniq
