# This is a Nomad Configuration file for server instances.
# I.e. the instances that assign Nomad jobs to run on client instances.

log_level = "WARN"

# Setup data dir
data_dir = "/tmp/nomad_server1"

# Enable the server
server {
  enabled = true

  # Self-elect, should be 3 or 5 for production
  bootstrap_expect = 3

  # Clean out old jobs quickly
  job_gc_threshold = "5m"

  # This is the IP address of the first server we provisioned
  retry_join = ["${nomad_lead_server_ip}:4648"]
}

consul {
  auto_advertise      = false
  server_auto_join    = false
  client_auto_join    = false
}

# Disabled because it slows Nomad down. We dream it will be fixed one day!
# telemetry {
#   publish_allocation_metrics = true
#   publish_node_metrics       = true
#   statsd_address = "${nomad_lead_server_ip}:8125"
# }
