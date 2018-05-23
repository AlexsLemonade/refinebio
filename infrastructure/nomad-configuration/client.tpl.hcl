# This is a Nomad Configuration file for client instances.
# I.e. the instances that actually run Nomad jobs.

log_level = "INFO"

# Setup data dir
data_dir = "/tmp/nomad_client1"

# Enable the client
client {
    enabled = true

    servers = ["${nomad_server_address}:4647"]
}
