#!/bin/bash

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

echo "Rebuilding Docker image while waiting for Nomad to come online."

# Start the nomad in development locally
nomad agent -bind $HOST_IP \
      -data-dir /home/kurt/Development/data_refinery/nomad_dir \
      -dev \
      > nomad.logs &

docker build -t wkurt/nomad-test:second -f workers/Dockerfile .

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

if [ ! -f workers/downloader.nomad -o ! -f workers/processor.nomad ]; then
    echo "One or more Job Specifications for Nomad are missing, building them from the template."
    ./workers/format_nomad_with_env.sh
fi

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 workers/downloader.nomad
nomad run -address http://$HOST_IP:4646 workers/processor.nomad
