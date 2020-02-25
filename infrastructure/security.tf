# The configuration in this file specifies AWS security groups and
# associated rules.

# This is the SSH key that can be used to ssh onto instances for
# debugging. Long term we may want to remove this entirely.
resource "aws_key_pair" "data_refinery" {
  key_name = "data-refinery-key-${var.user}-${var.stage}"
  public_key = "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDG1WFcvLLyLK7jZfhYnNBfbVwm259du2ig4tYcicA1d8d8I43LcWg2WYpd7EfNFH8LnJMDg632NcnQ0qzrpUG4zLGTcYufXm1Fm97J285iabzlUxfgSpbk5Ee1ioNCmqtPxEgy5lrt2xw0p3Rnbn0NvSKzwGU82/k/NCbxeKbaRpHLjz9TTcAdcZLugV7Syr8W+zWBqlCIMyC4ce4t8s/ecGbyacmRPdPqC9jUBC0guLHeQmlinINJIr+wMihxJ0B5Zcyokf4wXlQBPPcB89oO9L81nlApY6aK5JJrhkSN8M5+YOkdk6Xi4SZuJD5SLWbilKGPiCNiLPAnPw7m7Ual"
}

##
# Workers
##

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

# Allow the Nomad HTTP API to be accessible by the Foreman security group. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_nomad_from_foreman" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
  source_security_group_id = "${aws_security_group.data_refinery_foreman.id}"
}

# Allow the Nomad HTTP API to be accessible by the API security group. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_api_nomad" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
  source_security_group_id = "${aws_security_group.data_refinery_api.id}"
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

# Allow statsd updates from clients to lead server. See:
resource "aws_security_group_rule" "data_refinery_statsd" {
  type = "ingress"
  from_port = 8125
  to_port = 8125
  protocol = "udp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

# Allow SSH connections to Nomad instances. This is pretty much the
# only way to debug issues until we get the Nomad UI up and running.
# XXX: THIS DEFINITELY NEEDS TO BE REMOVED LONG TERM!!!!!!!!!!
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

##
# Database
##

resource "aws_security_group" "data_refinery_db" {
  name = "data-refinery_db-${var.user}-${var.stage}"
  description = "data_refinery_db-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-db-${var.user}-${var.stage}"
  }
}

resource "aws_security_group_rule" "data_refinery_db_outbound" {
  type = "egress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_security_group_rule" "data_refinery_db_pg_tcp" {
  type = "ingress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  source_security_group_id = "${aws_security_group.data_refinery_pg.id}"
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

##
# PGBouncer
##

resource "aws_security_group" "data_refinery_pg" {
  name = "data-refinery-pg-${var.user}-${var.stage}"
  description = "data_refinery_pg-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-pg-${var.user}-${var.stage}"
  }
}

resource "aws_security_group_rule" "data_refinery_pg_ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}


resource "aws_security_group_rule" "data_refinery_pg_workers_tcp" {
  type = "ingress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  source_security_group_id = "${aws_security_group.data_refinery_worker.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_foreman_tcp" {
  type = "ingress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  source_security_group_id = "${aws_security_group.data_refinery_foreman.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_api_tcp" {
  type = "ingress"
  from_port = 0
  to_port = 65535
  protocol = "tcp"
  source_security_group_id = "${aws_security_group.data_refinery_api.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_workers_icmp" {
  type = "ingress"
  from_port = -1
  to_port = -1
  protocol = "icmp"
  source_security_group_id = "${aws_security_group.data_refinery_worker.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_api_icmp" {
  type = "ingress"
  from_port = -1
  to_port = -1
  protocol = "icmp"
  source_security_group_id = "${aws_security_group.data_refinery_api.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_foreman_icmp" {
  type = "ingress"
  from_port = -1
  to_port = -1
  protocol = "icmp"
  source_security_group_id = "${aws_security_group.data_refinery_foreman.id}"
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

resource "aws_security_group_rule" "data_refinery_pg_outbound" {
  type = "egress"
  from_port = 0
  to_port = 65535
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}

##
# ElasticSearch
##

resource "aws_security_group" "data_refinery_es" {
  name = "data-refinery-es-${var.user}-${var.stage}"
  description = "data_refinery_es-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-es-${var.user}-${var.stage}"
  }

  # Wide open, but inside inside the VPC
  ingress {
      from_port = 0
      to_port = 0
      protocol = "-1"
      cidr_blocks = [ "0.0.0.0/0" ]
  }
}

##
# API
##

resource "aws_security_group" "data_refinery_api" {
  name = "data-refinery-api-${var.user}-${var.stage}"
  description = "data-refinery-api-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-api-${var.user}-${var.stage}"
  }
}

# XXX: THIS DEFINITELY NEEDS TO BE REMOVED LONG TERM!!!!!!!!!!
resource "aws_security_group_rule" "data_refinery_api_ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_api.id}"
}

resource "aws_security_group_rule" "data_refinery_api_http" {
  type = "ingress"
  from_port = 80
  to_port = 80
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_api.id}"
}

resource "aws_security_group_rule" "data_refinery_api_https" {
  type = "ingress"
  from_port = 443
  to_port = 443
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_api.id}"
}

resource "aws_security_group_rule" "data_refinery_api_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = "${aws_security_group.data_refinery_api.id}"
}

##
# Foreman
##

resource "aws_security_group" "data_refinery_foreman" {
  name = "data-refinery-foreman-${var.user}-${var.stage}"
  description = "data-refinery-foreman-${var.user}-${var.stage}"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-foreman-${var.user}-${var.stage}"
  }
}

# Allow the Nomad HTTP API to be accessible by this security group. See:
# https://www.nomadproject.io/guides/cluster/requirements.html#ports-used
resource "aws_security_group_rule" "data_refinery_worker_foreman" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_foreman.id}"
}

# XXX: THIS DEFINITELY NEEDS TO BE REMOVED LONG TERM!!!!!!!!!!
resource "aws_security_group_rule" "data_refinery_foreman_ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_foreman.id}"
}

# Necessary to retrieve
# https://s3.amazonaws.com/aws-cloudwatch/downloads/latest/awslogs-agent-setup.py
resource "aws_security_group_rule" "data_refinery_foreman_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = "${aws_security_group.data_refinery_foreman.id}"
}
