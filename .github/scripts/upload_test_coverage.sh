#!/bin/bash

# Script to upload code coverage to circleci

project=$1
if [[ $project == "" ]]
then
    echo "No project specified"
    exit 1
fi

if [[ $project == "workers" ]]
then
    # the workers project uses it's own test_volume directory
    test_volume="workers/test_volume"
else
    test_volume="test_volume"
fi

coverage_file="${test_volume}/coverage.xml"

if [[ ! -f $coverage_file ]]
then
    echo "Coverage file wasn't found, were the tests run before?"
    exit 0 # exit this script but don't fail the tests for this.
fi

output_file="${test_volume}/${project}_coverage.xml"

# In the test coverage report, all file paths are relative to each project
# folder. We need to be relative to the repo's root directory. That's why we
# append the project folder name to each file path in coverage.xml
sed "s/filename=\"/filename=\"$project\//g" "$coverage_file" > "$output_file"

# codecov.sh is located at https://codecov.io/bash
# we downloaded it for convenience
./.github/scripts/codecov.sh -f "$output_file" -Z -X s3 -F "$project"

rm -f "$coverage_file" "$output_file"
