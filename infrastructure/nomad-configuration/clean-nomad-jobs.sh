#!/bin/sh

TOKEN="$(curl -s -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")"

INSTANCE_ID="$(wget -q -O - http://169.254.169.254/latest/meta-data/instance-id)"
is_this_instance() {
    while read -r line; do
        job=$(nomad job inspect "$line" | jq '.Job.TaskGroups[].Tasks[] | select(.Constraints != null) | .Constraints[] | select(.LTarget == "${meta.volume_index}" and .RTarget == "'"$INSTANCE_ID"'")')
        [ -n "$job" ] && echo "$line"
    done
}

while sleep 5; do
    HTTP_CODE=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" -s -w '%{http_code}' -o /dev/null http://169.254.169.254/latest/meta-data/spot/instance-action)

    if [ "$HTTP_CODE" -eq 401 ] ; then
        TOKEN="$(curl -s -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 30")"
    elif [ "$HTTP_CODE" -eq 200 ] ; then
        nomad job status | grep parameterized | awk '{ print $1 }' | is_this_instance | xargs -L1 nomad stop
    else
        :
    fi

done
