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

# `coverage report -m` will always have an exit code of 0 which makes
# it seem like the test is passing. Therefore we store the exit code
# of running the tests as $exit_code, then report the coverage, and
# then exit with the appropriate code.
# This is done a function so arguments to the tests can be passed through.
run_tests_with_coverage () {
    echo "coverage run --source=\".\" manage.py test --no-input $@; exit_code=\$?; coverage report -m; exit \$exit_code"
}
