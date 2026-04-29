#!/bin/sh

# These are lists of docker images that we use. The actual names end
# up being <DOCKERHUB_REPO>/dr_<IMAGE_NAME> but this is useful for scripting.

# Sometimes we only want to work with the worker images.
export WORKER_IMAGES="smasher compendia illumina affymetrix salmon transcriptome no_op downloaders"

export ALL_IMAGES="base api_base api foreman $WORKER_IMAGES"

# This function checks whether a given docker image name ($1:$2)
# exists in Docker Hub or not using Docker Hub API V2. Based on:
# https://stackoverflow.com/questions/32113330/check-if-imagetag-combination-already-exists-on-docker-hub
docker_image_exists() {
    TOKEN=$(curl -s -H "Content-Type: application/json" -X POST \
        -d '{"username": "'"${DOCKER_USERNAME}"'", "password": "'"${DOCKER_PASSWORD}"'"}' \
        https://hub.docker.com/v2/users/login/ | jq -r .token)
    EXISTS=$(curl -s -H "Authorization: JWT ${TOKEN}" \
        "https://hub.docker.com/v2/repositories/$1/tags/?page_size=10000" |
        jq -r "[.results | .[] | .name == \"$2\"] | any" 2>/dev/null)
    test -n "$EXISTS" -a "$EXISTS" = true
}

get_branch_hash() {
    git rev-parse --abbrev-ref HEAD | shasum | awk '{print $1}'
}

get_docker_db_ip_address() {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' drdb 2>/dev/null
}

get_docker_es_ip_address() {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' dres 2>/dev/null
}

# A tag is linked to a commit hash, not a branch. A single commit hash
# can end up on multiple branches. So we first check to see if we're
# on master, then on dev, then error out because we should only deploy master or dev.
get_deploy_branch() {
    # Takes the version that is being deployed as its only parameter
    version="$1"

    if [ -z "$version" ]; then
        echo "You must pass the version to get_deploy_branch()."
    else
        is_master=$(git log origin/master --decorate=full | grep "$version" || true)
        is_dev=$(git log origin/dev --decorate=full | grep "$version" || true)

        # All dev versions should end with '-dev' or '-dev-hotfix' and all master versions should not.
        if [ -n "$is_master" ] && ! echo "$version" | grep -Eq "\-dev(\-hotfix)?$"; then
            echo "master"
        elif [ -n "$is_dev" ]; then
            echo "dev"
        else
            echo "unknown"
        fi
    fi
}

# Maps a branch ('master' or 'dev') to the corresponding deploy environment Dockerhub repo.
# Errors out for any other branch since deploys only happen from master and dev.
get_deploy_repo() {
    branch="$1"
    case "$branch" in
        master) echo "ccdl" ;;
        dev) echo "ccdlstaging" ;;
        *)
            echo "Error: get_deploy_repo requires 'master' or 'dev', got: '$branch'" >&2
            return 1
            ;;
    esac
}

# `coverage report -m` will always have an exit code of 0 which makes
# it seem like the test is passing. Therefore we store the exit code
# of running the tests as $exit_code, then report the coverage, and
# then exit with the appropriate code.
# This is done a function so arguments to the tests can be passed through.
run_tests_with_coverage() {
    COVERAGE="coverage run --source=\".\" manage.py test --settings=tests.settings --no-input $*; exit_code=\$?;"
    SAVE_REPORT="coverage xml -o data_store/coverage.xml;"
    PRINT_REPORT="coverage report -m;"
    RETURN="exit \$exit_code"

    echo "$COVERAGE $PRINT_REPORT $SAVE_REPORT $RETURN"
}

# Default Docker registry.
if [ -z "$DOCKERHUB_REPO" ]; then
    DOCKERHUB_REPO="ccdlstaging"
    export DOCKERHUB_REPO
fi

# Defaults to commit hash value for if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="$(get_branch_hash)"
    export SYSTEM_VERSION
fi
