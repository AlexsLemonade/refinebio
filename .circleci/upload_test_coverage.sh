#!/bin/sh

# Script to upload code coverage to circleci

project=$1
if [[ $project == "" ]]
then
    echo "No project specified"
    exit 1
fi

output="coverage_reports/${project}_coverage.xml"

sed "s/filename=\"/filename=\"$project\//g" coverage_reports/coverage.xml > $output

curl -s https://codecov.io/bash | bash -s -- -f "$output" -F unittests
