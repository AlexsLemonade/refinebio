#!/bin/bash

# Script to upload code coverage to circleci

project=$1
if [[ $project == "" ]]
then
    echo "No project specified"
    exit 1
fi

if [[ ! -f "test_volume/coverage.xml" ]]
then
    echo "Coverage file wasn't found, were the tests run before?"
    exit 1
fi

output="test_volume/${project}_coverage.xml"

# In the test coverage report, all file paths are absolute to each project
# folder. We need to be relative to refinebio/ that's why we append the project
# folder name to each file path in coverage.xml
sed "s/filename=\"/filename=\"$project\//g" test_volume/coverage.xml > $output

curl -s https://codecov.io/bash | bash -s -- -f "$output" -Z -F $project

rm -f test_volume/coverage.xml
rm -f $output
