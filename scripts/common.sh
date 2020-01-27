#!/bin/sh

# These are lists of docker images that we use. The actual names end
# up being <DOCKERHUB_REPO>/dr_<IMAGE_NAME> but this is useful for scripting.
export ALL_CCDL_IMAGES="smasher compendia illumina affymetrix salmon transcriptome no_op downloaders foreman api"
# Sometimes we only want to work with the worker images.
export CCDL_WORKER_IMAGES="smasher compendia illumina affymetrix salmon transcriptome no_op downloaders"

get_ip_address () {
    if [ "$(uname)" = "Linux" ]; then
        # If we're not connected to the internet return an error code.
        ip_out=$(ip route get 8.8.8.8)
        if [ $? -ne 0 ]; then
            return 1
        fi
        echo $ip_out | grep -oE 'src ([0-9]{1,3}\.){3}[0-9]{1,3}' | awk '{print $2; exit}'
    elif [ "$(uname)" = 'Darwin' ]; then # MacOS
         ip_out=$(ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2 | tail -1)
         if [ -z "$ip_out" ]; then
             return 1
         fi
         echo $ip_out

    fi
}

get_docker_db_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' drdb 2> /dev/null
}

get_docker_es_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' dres 2> /dev/null
}

# `coverage report -m` will always have an exit code of 0 which makes
# it seem like the test is passing. Therefore we store the exit code
# of running the tests as $exit_code, then report the coverage, and
# then exit with the appropriate code.
# This is done a function so arguments to the tests can be passed through.
run_tests_with_coverage () {
    COVERAGE="coverage run --source=\".\" manage.py test --no-input $*; exit_code=\$?;"
    SAVE_REPORT="coverage xml -o data_store/coverage.xml; sleep 1;"
    PRINT_REPORT="coverage report -m;"
    RETURN="exit \$exit_code"

    echo "$COVERAGE $SAVE_REPORT $PRINT_REPORT $RETURN"
}

# This function checks whether a given docker image name ($1:$CIRCLE_TAG)
# exists in Docker Hub or not using Docker Hub API V2. Based on:
# https://stackoverflow.com/questions/32113330/check-if-imagetag-combination-already-exists-on-docker-hub
docker_img_exists() {
    TOKEN=$(curl -s -H "Content-Type: application/json" -X POST \
                 -d '{"username": "'"${DOCKER_ID}"'", "password": "'"${DOCKER_PASSWD}"'"}' \
                 https://hub.docker.com/v2/users/login/ | jq -r .token)
    EXISTS=$(curl -s -H "Authorization: JWT ${TOKEN}" \
                  "https://hub.docker.com/v2/repositories/$1/tags/?page_size=10000" \
                 | jq -r "[.results | .[] | .name == \"$2\"] | any" 2> /dev/null)
    test -n "$EXISTS" -a "$EXISTS" = true
}

# A tag is linked to a commit hash, not a branch. A single commit hash
# can end up on multiple branches. So we first check to see if we're
# on master, then on dev, then error out because we should only deploy master or dev.
get_master_or_dev() {
    # Takes the version that is being deployed as its only parameter
    version="$1"

    if [ -z "$version" ]; then
        echo "You must pass the version to get_master_or_dev."
    else
        master_check=$(git log origin/master --decorate=full | grep "$version" || true)
        dev_check=$(git log origin/dev --decorate=full | grep "$version" || true)

        # All dev versions should end with '-dev' or '-dev-hotfix' and all master versions should not.
        if [ -n "$master_check" ] && ! echo "$version" | grep -Eq "\-dev(\-hotfix)?$"; then
            echo "master"
        elif [ -n "$dev_check" ] ; then
            echo "dev"
        else
            echo "unknown"
        fi
    fi
}

# Convenience function to export the NOMAD_ADDR environment variable,
# set to the address used for local development.
set_nomad_address() {
    NOMAD_ADDR="http://$(get_ip_address):4646"
    export NOMAD_ADDR
}

# Convenience function to export the NOMAD_ADDR environment variable,
# set to the address used for tests.
set_nomad_test_address() {
    NOMAD_ADDR="http://$(get_ip_address):5646"
    export NOMAD_ADDR
}
