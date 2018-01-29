#!/bin/bash

get_ip_address () {
    if [ `uname` == "Linux" ]; then
        echo $(ip route get 8.8.8.8 | awk '{print $NF; exit}')
    elif [ `uname` == 'Darwin' ]; then # MacOS
        echo $(ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2)
    fi
}

get_docker_db_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' drdb
}

get_docker_nomad_ip_address () {
    docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' nomad
}
