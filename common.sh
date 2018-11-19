#!/bin/bash

# These are lists of docker images that we use. The actual names end
# up being <DOCKERHUB_REPO>/dr_<IMAGE_NAME> but this is useful for scripting.
export ALL_CCDL_IMAGES="smasher illumina affymetrix salmon transcriptome no_op downloaders foreman api"
# Sometimes we only want to work with the worker images.
export CCDL_WORKER_IMAGES="smasher illumina affymetrix salmon transcriptome no_op downloaders"

get_ip_address () {
    if [ `uname` == "Linux" ]; then
        echo $(ip route get 8.8.8.8 | grep -oE 'src ([0-9]{1,3}\.){3}[0-9]{1,3}' | awk '{print $2; exit}')
    elif [ `uname` == 'Darwin' ]; then # MacOS
        echo $(ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2 | tail -1)
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

# A tag is linked to a commit hash, not a branch. A single commit hash
# can end up on multiple branches. So we first check to see if we're
# on master, then on dev, then error out because we should only deploy master or dev.
get_master_or_dev() {
    # Takes the version that is being deployed as its only parameter
    version=$1

    if [[ -z $version ]]; then
        echo "You must pass the version to get_master_or_dev."
        exit
    fi
    master_check=$(git log --decorate=full | grep -e $version -e origin/master || true)
    dev_check=$(git log --decorate=full | grep -e $version -e origin/dev || true)

    if [[ ! -z $master_check ]]; then
        echo "master"
    elif [[ ! -z $dev_check ]]; then
        echo "dev"
    else
        echo "Why in the world was update_docker_img.sh called from a branch other than dev or master?!?!?"
        exit
    fi
}
