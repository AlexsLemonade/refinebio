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
        -d '{"username": "'"${DOCKER_ID}"'", "password": "'"${DOCKER_PASSWD}"'"}' \
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

# Create docker-container/desktop-linux Docker builder if none provided.
# Set the builder as currently used.
set_up_docker_builder() {
    if [ -z "$DOCKER_BUILDER" ]; then
        DOCKER_BUILDER="refinebio_local_builder"
        echo "Creating Docker builder $DOCKER_BUILDER."
        docker buildx create \
            --driver=docker-container \
            --name="$DOCKER_BUILDER" \
            --platform=linux/amd64 2>/dev/null ||
            true
    fi

    docker buildx use "$DOCKER_BUILDER"
    echo "Using Docker builder $DOCKER_BUILDER:"
    docker buildx inspect
}

update_docker_image() {
    DOCKERHUB_REPO="$1"
    IMAGE_NAME="$2"
    SYSTEM_VERSION="$3"
    DOCKER_FILE_PATH="$4"
    DOCKER_ACTION="$5"

    if [ -z "$DOCKER_ACTION" ]; then
        DOCKER_ACTION="--load"
    fi

    CCDL_STAGING_IMAGE="ccdlstaging/dr_$IMAGE_NAME"
    DOCKERHUB_IMAGE="$DOCKERHUB_REPO/dr_$IMAGE_NAME"
    CACHE_FROM_CCDL_STAGING_LATEST="type=registry,ref=${CCDL_STAGING_IMAGE}_cache:latest"
    CACHE_FROM_CCDL_STAGING_VERSION="type=registry,ref=${CCDL_STAGING_IMAGE}_cache:$SYSTEM_VERSION"
    CACHE_FROM_LATEST="type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest"
    CACHE_FROM_VERSION="type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION"
    CACHE_TO_LATEST="type=registry,ref=${DOCKERHUB_IMAGE}_cache:latest,mode=max"
    CACHE_TO_VERSION="type=registry,ref=${DOCKERHUB_IMAGE}_cache:$SYSTEM_VERSION,mode=max"


    set_up_docker_builder

    docker buildx build \
        --build-arg DOCKERHUB_REPO="$DOCKERHUB_REPO" \
        --build-arg SYSTEM_VERSION="$SYSTEM_VERSION" \
        --cache-from "$CACHE_FROM_CCDL_STAGING_LATEST" \
        --cache-from "$CACHE_FROM_CCDL_STAGING_VERSION" \
        --cache-from "$CACHE_FROM_LATEST" \
        --cache-from "$CACHE_FROM_VERSION" \
        --cache-to "$CACHE_TO_LATEST" \
        --cache-to "$CACHE_TO_VERSION" \
        --file "$DOCKER_FILE_PATH" \
        --platform linux/amd64 \
        --tag "$DOCKERHUB_IMAGE:latest" \
        --tag "$DOCKERHUB_IMAGE:$SYSTEM_VERSION" \
        "$DOCKER_ACTION" \
        .
}

# Default Docker registry.
if [ -z "$DOCKERHUB_REPO" ]; then
    export DOCKERHUB_REPO="ccdlstaging"
fi

# Defaults to commit hash value for if we're not running in the cloud.
if [ -z "$SYSTEM_VERSION" ]; then
    SYSTEM_VERSION="$(get_branch_hash)"
fi
