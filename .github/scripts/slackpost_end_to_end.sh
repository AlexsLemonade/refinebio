#!/bin/bash

if [[ $ENGAGEMENT_WEBHOOK == "" ]]
then
    echo "No webhook url. Set ENGAGEMENT_WEBHOOK in the environment variables if you want to be notified of deploys on slack"
    exit 0
fi

# ------------
channel=$1
if [[ $channel == "" ]]
then
        echo "No channel specified"
        exit 1
fi

# ------------
shift
username=$1
if [[ $username == "" ]]
then
        echo "No username specified"
        exit 1
fi

text="The end-to-end tests passed in the staging stack!!!"

escapedText=$(echo "$text" | sed 's/"/\"/g' | sed "s/'/\'/g" )

json="{\"channel\": \"$channel\", \"username\":\"$username\", \"icon_emoji\":\":tada:\", \"attachments\":[{\"color\":\"danger\" , \"text\": \"$escapedText\"}]}"

curl -s -d "payload=$json" "$ENGAGEMENT_WEBHOOK"
