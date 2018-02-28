# The configuration in this file specifies AWS security groups and
# associated rules.

# This is the SSH key that can be used to ssh onto instances for
# debugging. Long term we may want to remove this entirely.
resource "aws_key_pair" "data_refinery" {
  key_name = "data-refinery-key-${var.user}-${var.stage}"
  public_key = "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDG1WFcvLLyLK7jZfhYnNBfbVwm259du2ig4tYcicA1d8d8I43LcWg2WYpd7EfNFH8LnJMDg632NcnQ0qzrpUG4zLGTcYufXm1Fm97J285iabzlUxfgSpbk5Ee1ioNCmqtPxEgy5lrt2xw0p3Rnbn0NvSKzwGU82/k/NCbxeKbaRpHLjz9TTcAdcZLugV7Syr8W+zWBqlCIMyC4ce4t8s/ecGbyacmRPdPqC9jUBC0guLHeQmlinINJIr+wMihxJ0B5Zcyokf4wXlQBPPcB89oO9L81nlApY6aK5JJrhkSN8M5+YOkdk6Xi4SZuJD5SLWbilKGPiCNiLPAnPw7m7Ual"
}

# This is a security group for Data Refinery Workers, which currently
# includes the Nomad Server nodes as well.
resource "aws_security_group" "data_refinery_worker" {
  name = "data-refinery-worker-${var.user}-${var.stage}"
  description = "data-refinery-worker-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-worker-${var.user}-${var.stage}"
  }
}

# Allow HTTP connections from this security group to itself.
resource "aws_security_group_rule" "data_refinery_worker_custom" {
  type = "ingress"
  from_port = 8000
  to_port = 8000
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow the Nomad HTTP API to be accessible by this security group. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_worker_nomad" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow Nomad RPC calls to go through to each other. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_worker_nomad_rpc" {
  type = "ingress"
  from_port = 4647
  to_port = 4647
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow Nomad Servers to gossip with each other over TCP. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_worker_nomad_serf_tcp" {
  type = "ingress"
  from_port = 4648
  to_port = 4648
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow Nomad Servers to gossip with each other over UDP. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_worker_nomad_serf_udp" {
  type = "ingress"
  from_port = 4648
  to_port = 4648
  protocol = "udp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow SSH connections to Nomad instances. This is pretty much the
# only way to debug issues until we get the Nomad UI up and running.
# THIS DEFINITELY NEEDS TO BE REMOVED LONG TERM!!!!!!!!!!
resource "aws_security_group_rule" "data_refinery_worker_ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow outbound requests from Nomad so they can actually do useful
# things. Long term we will probably only want to allow these to come
# from the Nomad Client instances, but we do not yet have separate
# security groups for servers and clients.
resource "aws_security_group_rule" "data_refinery_worker_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_security_group" "data_refinery_db" {
  name = "data-refinery_db-${var.user}-${var.stage}"
  description = "data_refinery_db-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-db-${var.user}-${var.stage}"
  }
}

resource "aws_security_group_rule" "data_refinery_db_postgres_lab" {
  type = "ingress"
  from_port = 5432
  to_port = 5432
  protocol = "tcp"
  cidr_blocks = ["165.123.67.153/32", "165.123.67.173/32"] # Lab machines
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_security_group_rule" "data_refinery_db_workers_tcp" {
  type = "ingress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  source_security_group_id = "${aws_security_group.data_refinery_worker.id}"
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_security_group_rule" "data_refinery_db_workers_icmp" {
  type = "ingress"
  from_port = -1
  to_port = -1
  protocol = "icmp"
  source_security_group_id = "${aws_security_group.data_refinery_worker.id}"
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_security_group_rule" "data_refinery_db_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}
