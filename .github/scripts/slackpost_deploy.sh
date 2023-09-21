#!/bin/bash

if [[ $SLACK_WEBHOOK_URL == "" ]]; then
    echo "No webhook url. Set SLACK_WEBHOOK_URL in the environment variables if you want to be notified of deploys on slack"
    exit 0
fi

# ------------
channel=$1
if [[ $channel == "" ]]; then
    echo "No channel specified"
    exit 1
fi

# ------------
shift
username=$1
if [[ $username == "" ]]; then
    echo "No username specified"
    exit 1
fi

# ------------
master_check=$(git branch --contains "tags/$CI_TAG" | grep '^  master$' || true)
dev_check=$(git branch --contains "tags/$CI_TAG" | grep '^  dev$' || true)

if [[ -n $master_check ]]; then
    CI_BRANCH=master
elif [[ -n $dev_check ]]; then
    CI_BRANCH=dev
fi

text="New deployment! Woo! $CI_USERNAME: $CI_BRANCH $CI_TAG"

escapedText=$(echo "$text" | sed 's/"/\"/g' | sed "s/'/\'/g")

json="{\"channel\": \"$channel\", \"username\":\"$username\", \"icon_emoji\":\":tada:\", \"attachments\":[{\"color\":\"danger\" , \"text\": \"$escapedText\"}]}"

curl -s -d "payload=$json" "$SLACK_WEBHOOK_URL"
