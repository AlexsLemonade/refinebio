# This is a Nomad Configuration file for client instances.
# I.e. the instances that actually run Nomad jobs.

# Increase log verbosity for now.
log_level = "DEBUG"

# Setup data dir
data_dir = "/tmp/nomad_client1"

# Enable the client
client {
    enabled = true

    # This doesn't work because getting the IP of the Nomad server is non-trivial.
    # TODO: Have the client-instance-user-data.sh be templated by terraform such that it
    # writes the server's IP address to an environment variable so that it can be used here.
    servers = ["54.227.72.146:4647"]
}
