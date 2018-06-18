#!/bin/bash

get_ip_address () {
    if [ `uname` == "Linux" ]; then
        echo $(dig +short myip.opendns.com @resolver1.opendns.com)
    elif [ `uname` == 'Darwin' ]; then # MacOS
        echo $(ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2)
    fi
}

get_docker_db_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' drdb 2> /dev/null
}

# `coverage report -m` will always have an exit code of 0 which makes
# it seem like the test is passing. Therefore we store the exit code
# of running the tests as $exit_code, then report the coverage, and
# then exit with the appropriate code.
# This is done a function so arguments to the tests can be passed through.
run_tests_with_coverage () {
    echo "coverage run --source=\".\" manage.py test --no-input $@; exit_code=\$?; coverage report -m; exit \$exit_code"
}

# This function checks whether a given docker image name ($1:$CIRCLE_TAG)
# exists in Docker Hub or not using Docker Hub API V2. Based on:
# https://stackoverflow.com/questions/32113330/check-if-imagetag-combination-already-exists-on-docker-hub
function docker_img_exists() {
    TOKEN=$(curl -s -H "Content-Type: application/json" -X POST \
                 -d '{"username": "'${DOCKER_ID}'", "password": "'${DOCKER_PASSWD}'"}' \
                 https://hub.docker.com/v2/users/login/ | jq -r .token)
    EXISTS=$(curl -s -H "Authorization: JWT ${TOKEN}" \
                  https://hub.docker.com/v2/repositories/$1/tags/?page_size=10000 \
             | jq -r "[.results | .[] | .name == \"$2\"] | any" 2> /dev/null)
    test -n "$EXISTS" -a "$EXISTS" = true
}
