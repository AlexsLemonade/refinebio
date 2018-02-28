# The configuration contained in this file specifies the AWS instances
# we'll need. These include EC2 instances and an RDS instance.

data "aws_ami" "ubuntu" {
  most_recent = true

  filter {
    name   = "name"
    # This is a HVM, EBS backed SSD Ubuntu LTS AMI with Docker version 17.12 on it.
    values = ["ubuntu-16.04-docker-17.12.0-ce-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }
}

resource "aws_instance" "data_refinery_worker_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.xlarge"
  availability_zone = "us-east-1a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  user_data = "${file("server-instance-user-data.sh")}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "data-refinery-1-${var.user}-${var.stage}"
  }

  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }

  ebs_block_device = {
    device_name = "/dev/xvdcz"
    volume_type = "gp2"
    volume_size = 40
  }
}

resource "aws_instance" "data_refinery_worker_2" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.xlarge"
  availability_zone = "us-east-1b"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  depends_on = ["aws_internet_gateway.data_refinery", "aws_instance.data_refinery_worker_1"]
  user_data = "${file("client-instance-user-data.sh")}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "data-refinery-2-${var.user}-${var.stage}"
  }

  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }

  ebs_block_device = {
    device_name = "/dev/xvdcz"
    volume_type = "gp2"
    volume_size = 40
  }
}

variable "database_password" {}

resource "aws_db_instance" "postgres-db" {
  identifier = "data-refinery-${var.user}-${var.stage}"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.5.4"
  instance_class = "db.t2.micro"
  name = "data-refinery-${var.user}-${var.stage}"
  username = "data_refinery_user"
  password = "${var.database_password}"
  db_subnet_group_name = "${aws_db_subnet_group.data_refinery.name}"
  skip_final_snapshot = true
  vpc_security_group_ids = ["${aws_security_group.data_refinery_db.id}"]
  multi_az = true
  publicly_accessible = true
}
