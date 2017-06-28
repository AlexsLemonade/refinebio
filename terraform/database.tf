provider "aws" {
  region = "us-east-1"
}

variable "database_password" {}

resource "aws_vpc" "data_refinery_vpc" {
  cidr_block = "10.0.0.0/16"
  enable_dns_support = true
  enable_dns_hostnames = true

  tags {
    Name = "data-refinery"
  }
}

resource "aws_subnet" "data_refinery_1a" {
  availability_zone = "us-east-1a"
  cidr_block = "10.0.0.0/17"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"
  map_public_ip_on_launch = true

  tags {
    Name = "data-refinery-1a"
  }
}

resource "aws_subnet" "data_refinery_1b" {
  availability_zone = "us-east-1b"
  cidr_block = "10.0.128.0/17"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"
  # Unsure if this should be set to true
  map_public_ip_on_launch = true

  tags {
    Name = "data-refinery-1b"
  }
}

resource "aws_internet_gateway" "data_refinery" {
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags = {
    Name = "data-refinery"
  }
}

# Note: this is a insecure practice long term, however it's
# necessary to access it from lab machines.
resource "aws_route_table" "data_refinery" {
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = "${aws_internet_gateway.data_refinery.id}"
  }

  tags {
    Name = "data_refinery"
  }
}

resource "aws_route_table_association" "data_refinery_1a" {
  subnet_id      = "${aws_subnet.data_refinery_1a.id}"
  route_table_id = "${aws_route_table.data_refinery.id}"
}

resource "aws_route_table_association" "data_refinery_1b" {
  subnet_id      = "${aws_subnet.data_refinery_1b.id}"
  route_table_id = "${aws_route_table.data_refinery.id}"
}

resource "aws_db_subnet_group" "data_refinery" {
  name = "data_refinery"
  subnet_ids = ["${aws_subnet.data_refinery_1a.id}", "${aws_subnet.data_refinery_1b.id}"]

  tags {
    Name = "Data Refinery DB Subnet"
  }
}

resource "aws_security_group" "data_refinery_db" {
  name = "data_refinery_db"
  description = "data_refinery_db"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-db"
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

resource "aws_security_group_rule" "data_refinery_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_db_instance" "postgres-db" {
  identifier = "data-refinery"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.5.4"
  instance_class = "db.t2.large"
  name = "data_refinery"
  username = "data_refinery_user"
  password = "${var.database_password}"
  db_subnet_group_name = "${aws_db_subnet_group.data_refinery.name}"
  skip_final_snapshot = true
  vpc_security_group_ids = ["${aws_security_group.data_refinery_db.id}"]
  multi_az = true
  publicly_accessible = true
}

resource "aws_sqs_queue" "data_refinery_queue" {
  name = "data-refinery-queue"
  receive_wait_time_seconds = 10
}
