#!/bin/bash

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

# Start the nomad in development locally
nomad agent -bind $HOST_IP \
      -data-dir /home/kurt/Development/data_refinery/nomad_dir \
      -dev \
      > nomad.logs &

echo "Waiting for Nomad to come online."

check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  http://$HOST_IP:4646/v1/status/leader)
}

# Wait for Nomad to get started.
nomad_status=$(check_nomad_status)
while [ $nomad_status != "200" ]; do
    sleep 0.5
    nomad_status=$(check_nomad_status)
done

echo "Nomad is online. Registering jobs."

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 downloader.nomad
nomad run -address http://$HOST_IP:4646 processor.nomad
