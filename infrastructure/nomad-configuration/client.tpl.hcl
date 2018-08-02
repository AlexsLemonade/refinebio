# This is a Nomad Configuration file for client instances.
# I.e. the instances that actually run Nomad jobs.

log_level = "WARN"

# Setup data dir
data_dir = "/tmp/nomad_client1"

# Enable the client
client {
    enabled = true

    servers = ["${nomad_lead_server_ip}:4647"]
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
