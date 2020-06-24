#!/bin/bash

# Usage: slackpost "<channel>" "<username>" "<message>"

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

# ------------
master_check=$(git branch --contains "tags/$CI_TAG" | grep '^  master$' || true)
dev_check=$(git branch --contains "tags/$CI_TAG" | grep '^  dev$' || true)

if [[ ! -z $master_check ]]; then
    CI_BRANCH=master
elif [[ ! -z $dev_check ]]; then
    CI_BRANCH=dev
fi

text="New deployment! Woo! $CI_USERNAME: $CI_BRANCH $CI_TAG"

escapedText=$(echo "$text" | sed 's/"/\"/g' | sed "s/'/\'/g" )

json="{\"channel\": \"$channel\", \"username\":\"$username\", \"icon_emoji\":\":veerapan:\", \"attachments\":[{\"color\":\"danger\" , \"text\": \"$escapedText\"}]}"

curl -s -d "payload=$json" "$ENGAGEMENT_WEBHOOK"
