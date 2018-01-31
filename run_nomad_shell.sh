#! /bin/bash

source common.sh

docker run \
       --net data-refinery \
       --add-host=nomad:$(get_docker_nomad_ip_address) \
       --link nomad:nomad \
       miserlou/nomad $1 -address http://nomad:4646 "${@:2}"
