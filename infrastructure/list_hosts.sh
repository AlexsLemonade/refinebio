#! /bin/bash
# Use like:
# ./list_hosts.sh data-refinery-key-circleci-prod > hosts
# ./connect_cluster.sh

if [[ $# -eq 0 ]] ; then
    echo "Hey, you need to supply a user!"
    exit 0
fi
if [[ $# -eq 1 ]] ; then
    AWS_REGION=us-east-1
else
    AWS_REGION=$2
fi

aws ec2 describe-instances --region="$AWS_REGION" \
                           --filters "Name=tag:User,Values=$1" \ 
                           | grep PublicDnsName | tr -d '"' \
                           | sed "s/PublicDnsName: //g" \
                           | tr -d "," | awk '{$1=$1};1' | uniq
