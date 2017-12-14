#!/bin/bash

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

# Start the nomad in development locally
nomad agent -bind $HOST_IP -data-dir /home/kurt/Development/data_refinery/nomad_dir -dev > nomad.logs &

# Hacky way to wait for Nomad to get started
sleep 5

# Register the jobs for dispatching.
nomad run -address http://$HOST_IP:4646 downloader.nomad
nomad run -address http://$HOST_IP:4646 processor.nomad
