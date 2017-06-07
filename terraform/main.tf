provider "aws" {
  region = "us-east-1"
}

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

resource "aws_iam_instance_profile" "ecs_instance_profile" {
  name  = "data-refinery-ecs-instance-profile"
  roles = ["${aws_iam_role.ecs_instance.name}"]
}

resource "aws_iam_role" "ecs_instance" {
  name = "data-refinery-ecs-instance"

  assume_role_policy = <<EOF
{
  "Version": "2008-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF
}

resource "aws_iam_policy_attachment" "ecs" {
  name = "AmazonEC2ContainerServiceforEC2Role"
  roles = ["${aws_iam_role.ecs_instance.name}"]

  # The following can be found here:
  # https://console.aws.amazon.com/iam/home?region=us-east-1#/policies/arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_key_pair" "data_refinery" {
  key_name = "data-refinery-key"
  public_key = "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDG1WFcvLLyLK7jZfhYnNBfbVwm259du2ig4tYcicA1d8d8I43LcWg2WYpd7EfNFH8LnJMDg632NcnQ0qzrpUG4zLGTcYufXm1Fm97J285iabzlUxfgSpbk5Ee1ioNCmqtPxEgy5lrt2xw0p3Rnbn0NvSKzwGU82/k/NCbxeKbaRpHLjz9TTcAdcZLugV7Syr8W+zWBqlCIMyC4ce4t8s/ecGbyacmRPdPqC9jUBC0guLHeQmlinINJIr+wMihxJ0B5Zcyokf4wXlQBPPcB89oO9L81nlApY6aK5JJrhkSN8M5+YOkdk6Xi4SZuJD5SLWbilKGPiCNiLPAnPw7m7Ual"
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

resource "aws_security_group" "data_refinery_worker" {
  name = "data-refinery-worker"
  description = "data-refinery-worker"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags {
    Name = "data-refinery-worker"
  }
}

# resource "aws_security_group_rule" "cognoma-service-http" {
#   type = "ingress"
#   from_port = 80
#   to_port = 80
#   protocol = "tcp"
#   self = true
#   security_group_id = "${aws_security_group.cognoma-service.id}"
# }

resource "aws_security_group_rule" "data_refinery_worker_custom" {
  type = "ingress"
  from_port = 8000
  to_port = 8000
  protocol = "tcp"
  self = true
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_security_group_rule" "data_refinery_worker_ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_security_group_rule" "data_refinery_worker_outbound" {
  type = "egress"
  from_port = 0
  to_port = 0
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_instance" "data_refinery_worker_1" {
  ami = "ami-275ffe31"
  instance_type = "m4.large"
  availability_zone = "us-east-1a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  user_data = "${file("instance-user-data.sh")}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "data-refinery-1"
  }
}

resource "aws_instance" "data_refinery_worker_2" {
  ami = "ami-275ffe31"
  instance_type = "m4.large"
  availability_zone = "us-east-1b"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  user_data = "${file("instance-user-data.sh")}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "data-refinery-2"
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

resource "aws_elb" "data_refinery_worker" {
  name = "data-refinery-worker"
  subnets = [
    "${aws_subnet.data_refinery_1a.id}",
    "${aws_subnet.data_refinery_1b.id}"]

  listener {
    instance_port = 8000
    instance_protocol = "http"
    lb_port = 80
    lb_protocol = "http"
  }

  health_check {
    healthy_threshold = 10
    unhealthy_threshold = 2
    timeout = 5
    target = "TCP:8000"
    interval = 30
  }

  security_groups = ["${aws_security_group.data_refinery_worker.id}"]

  instances = [
    "${aws_instance.data_refinery_worker_1.id}",
    "${aws_instance.data_refinery_worker_2.id}",
  ]

  cross_zone_load_balancing = true
  idle_timeout = 60
  connection_draining = true
  connection_draining_timeout = 400
  internal = true
}
