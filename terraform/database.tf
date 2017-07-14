# The configuration contained in this file specifies AWS resources
# related to the database of The Data Refinery such as the RDS
# instance itself along with supporting resources such as subnets and
# security groups.

variable "database_password" {}

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

resource "aws_db_instance" "postgres-db" {
  identifier = "data-refinery"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.5.4"
  instance_class = "db.t2.micro"
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
