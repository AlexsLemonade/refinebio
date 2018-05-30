# This is a Nomad Configuration file for server instances.
# I.e. the instances that assign Nomad jobs to run on client instances.

log_level = "INFO"

# Setup data dir
data_dir = "/tmp/nomad_server1"

# Enable the server
server {
    enabled = true

    # Only 1 for the lead server, 3 for the others. Tx Kurt.
    bootstrap_expect = 1
}

consul {
  auto_advertise      = false
  server_auto_join    = false
  client_auto_join    = false
}

