# This is a Nomad Configuration file for server instances.
# I.e. the instances that assign Nomad jobs to run on client instances.

log_level = "INFO"

# Setup data dir
data_dir = "/tmp/nomad_server1"

# Enable the server
server {
  enabled = true

  # Self-elect, should be 3 or 5 for production
  bootstrap_expect = 3

  # This is the IP address of the first server we provisioned
  retry_join = ["${nomad_lead_server_ip}:4648"]
}

consul {
  auto_advertise      = false
  server_auto_join    = false
  client_auto_join    = false
}

telemetry {
  publish_allocation_metrics = true
  publish_node_metrics       = true
  statsd_address = "${nomad_lead_server_ip}:8125"
}
