#!/bin/bash

get_ip_address () {
    if [ `uname` == "Linux" ]; then
        echo $(ip route get 8.8.8.8 | awk '{print $NF; exit}')
    elif [ `uname` == 'Darwin' ]; then # MacOS
        echo $(ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2)
    fi
}

get_docker_db_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' drdb 2> /dev/null
}

get_docker_nomad_ip_address () {
    # If we aren't running Nomad in a container then it should be
    # running on the host, which means that the host's IP can be used
    # for Nomad.
    # The echo is in this command because without it an extra newline is output.
    temp=$(docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' nomad 2> /dev/null) \
        && echo $temp \
        || get_ip_address
}

get_nomad_link_option () {
    # If Nomad is running in a container then it won't have the same IP
    # address as the host and we'll need to link to it. If Nomad is
    # running on the host then we shouldn't try to link to the
    # non-existant container.
    if [ "$HOST_IP" != "$NOMAD_HOST_IP" ]; then
        echo "--link nomad:nomad"
    fi
}

# `coverage report -m` will always have an exit code of 0 which makes
# it seem like the test is passing. Therefore we store the exit code
# of running the tests as $exit_code, then report the coverage, and
# then exit with the appropriate code.
# This is done a function so arguments to the tests can be passed through.
run_tests_with_coverage () {
    echo "coverage run --source=\".\" manage.py test --no-input $@; exit_code=\$?; coverage report -m; exit \$exit_code"
}
