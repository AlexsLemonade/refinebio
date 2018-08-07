#!/bin/bash

# Usage: slackpost "<channel>" "<username>" "<message>"

# ------------
webhook_url="https://hooks.slack.com/services/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"

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

# ------------
master_check=$(git branch --contains tags/$CIRCLE_TAG | grep '^  master$')
dev_check=$(git branch --contains tags/$CIRCLE_TAG | grep '^  dev$')

if [[ ! -z $master_check ]]; then
    CIRCLE_BRANCH=master
elif [[ ! -z $dev_check ]]; then
    CIRCLE_BRANCH=dev
fi

text="New deployment! Woo! $CIRCLE_PROJECT_USERNAME: $CIRCLE_PULL_REQUEST $CIRCLE_BRANCH $CIRCLE_TAG"

escapedText=$(echo $text | sed 's/"/\"/g' | sed "s/'/\'/g" )

json="{\"channel\": \"$channel\", \"username\":\"$username\", \"icon_emoji\":\":veerapan:\", \"attachments\":[{\"color\":\"danger\" , \"text\": \"$escapedText\"}]}"

curl -s -d "payload=$json" "$webhook_url"
