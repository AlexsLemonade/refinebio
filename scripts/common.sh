#!/bin/sh

# These are lists of docker images that we use. The actual names end
# up being <DOCKERHUB_REPO>/dr_<IMAGE_NAME> but this is useful for scripting.
export ALL_IMAGES="smasher compendia illumina affymetrix salmon transcriptome no_op downloaders foreman api"
# Sometimes we only want to work with the worker images.
export WORKER_IMAGES="smasher compendia illumina affymetrix salmon transcriptome no_op downloaders"

# Default Docker registry.
if [ -z "$DOCKERHUB_REPO" ]; then
    export DOCKERHUB_REPO="ccdlstaging"
fi

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
        elif [ -n "$dev_check" ]; then
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
        echo "Setting up refine.bio Docker builder."

        if test "$GITHUB_ACTION"; then
            echo "$INSTANCE_SSH_KEY" >infrastructure/data-refinery-key.pem
            chmod 600 infrastructure/data-refinery-key.pem

            # shellcheck disable=SC2046
            eval $(ssh-agent)
            ssh-add infrastructure/data-refinery-key.pem

            if [ ! -d ~/.ssh ]; then
                mkdir -m 700 ~/.ssh
            fi
            cat >~/.ssh/config <<EOF
Host "${DEPLOY_IP_ADDRESS}"
    StrictHostKeyChecking no
    UserKnownHostsFile=/dev/null
EOF
            DOCKER_BUILDER="refinebio_remote_builder"
            echo "Creating Docker builder $DOCKER_BUILDER."
            docker buildx create \
                --driver=docker-container \
                --name="$DOCKER_BUILDER" \
                --platform=linux/amd64 \
                "ssh://ubuntu@${DEPLOY_IP_ADDRESS}" 2>/dev/null ||
                true
        else
            DOCKER_BUILDER="refinebio_local_builder"
            echo "Creating Docker builder $DOCKER_BUILDER."
            docker buildx create \
                --driver=docker-container \
                --name="$DOCKER_BUILDER" \
                --platform=linux/amd64 2>/dev/null ||
                true
        fi
    fi

    echo "Using Docker builder $DOCKER_BUILDER."
    docker buildx use "$DOCKER_BUILDER"
}
