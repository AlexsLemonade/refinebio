# The configuration contained in this file specifies high level AWS
# resources such as VPC, subnets, and security groups needed for The
# Data Refinery. It also is the current home of the EC2 instances and
# the IAM permissions/roles for them.

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
  role = "${aws_iam_role.ecs_instance.name}"
}

resource "aws_iam_role" "ecs_instance" {
  name = "data-refinery-ecs-instance"

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html
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

resource "aws_iam_role_policy_attachment" "ecs" {
  role = "${aws_iam_role.ecs_instance.name}"

  # The following can be found here:
  # https://console.aws.amazon.com/iam/home?region=us-east-1#/policies/arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_policy" "s3_access_policy" {
  name = "data-refinery-s3-access-policy"
  description = "Allows S3 Permissions."

  # Policy text based off of:
  # http://docs.aws.amazon.com/AmazonS3/latest/dev/example-bucket-policies.html
  policy = <<EOF
{
   "Version":"2012-10-17",
   "Statement":[
      {
         "Effect":"Allow",
         "Action":[
            "s3:ListAllMyBuckets"
         ],
         "Resource":"arn:aws:s3:::*"
      },
      {
         "Effect":"Allow",
         "Action":[
            "s3:ListBucket",
            "s3:GetBucketLocation"
         ],
         "Resource":"arn:aws:s3:::data-refinery"
      },
      {
         "Effect":"Allow",
         "Action":[
            "s3:PutObject",
            "s3:GetObject",
            "s3:DeleteObject"
         ],
         "Resource":"arn:aws:s3:::data-refinery/*"
      }
   ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "s3" {
  role = "${aws_iam_role.ecs_instance.name}"
  policy_arn = "${aws_iam_policy.s3_access_policy.arn}"
}

resource "aws_iam_policy" "cloudwatch_policy" {
  name = "data-refinery-cloudwatch-policy"
  description = "Allows Cloudwatch Permissions."


  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/iam-identity-based-access-control-cwl.html
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "logs:CreateLogGroup",
                "logs:CreateLogStream",
                "logs:PutLogEvents",
                "logs:DescribeLogStreams"
            ],
            "Resource": [
                "arn:aws:logs:*:*:*"
            ]
        }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "cloudwatch" {
  role = "${aws_iam_role.ecs_instance.name}"
  policy_arn = "${aws_iam_policy.cloudwatch_policy.arn}"
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
  ami = "ami-04351e12"
  instance_type = "t2.xlarge"
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

  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }

  ebs_block_device = {
    device_name = "/dev/xvdcz"
    volume_type = "gp2"
    volume_size = 40
  }

  provisioner "remote-exec" {
    # Commands copied from
    # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/using_cloudwatch_logs.html
    inline = [
      "sudo yum install -y awslogs",
      "sudo yum install -y jq",
    ]

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
  }

  provisioner "file" {
    source = "conf/awslogs.conf"
    destination = "/home/ec2-user/awslogs.conf"

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
  }

  provisioner "remote-exec" {
    # Commands copied from
    # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/using_cloudwatch_logs.html
    inline = [
      "sudo mv /home/ec2-user/awslogs.conf /etc/awslogs/awslogs.conf",
      "cluster=$(curl -s http://localhost:51678/v1/metadata | jq -r '. | .Cluster')",
      "sudo sed -i -e \"s/{cluster}/$cluster/g\" /etc/awslogs/awslogs.conf",
      "container_instance_id=$(curl -s http://localhost:51678/v1/metadata | jq -r '. | .ContainerInstanceArn' | awk -F/ '{print $2}' )",
      "sudo sed -i -e \"s/{container_instance_id}/$container_instance_id/g\" /etc/awslogs/awslogs.conf",
      "sudo service awslogs start",
      "sudo chkconfig awslogs on",
    ]

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
  }
}

resource "aws_instance" "data_refinery_worker_2" {
  ami = "ami-04351e12"
  instance_type = "t2.xlarge"
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

  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }

  ebs_block_device = {
    device_name = "/dev/xvdcz"
    volume_type = "gp2"
    volume_size = 40
  }

  provisioner "remote-exec" {
    # Commands copied from
    # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/using_cloudwatch_logs.html
    inline = [
      "sudo yum install -y awslogs",
      "sudo yum install -y jq",
    ]

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
  }

  provisioner "file" {
    source = "conf/awslogs.conf"
    destination = "/home/ec2-user/awslogs.conf"

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
  }

  provisioner "remote-exec" {
    # Commands copied from
    # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/using_cloudwatch_logs.html
    inline = [
      "sudo mv /home/ec2-user/awslogs.conf /etc/awslogs/awslogs.conf",
      "cluster=$(curl -s http://localhost:51678/v1/metadata | jq -r '. | .Cluster')",
      "sudo sed -i -e \"s/{cluster}/$cluster/g\" /etc/awslogs/awslogs.conf",
      "container_instance_id=$(curl -s http://localhost:51678/v1/metadata | jq -r '. | .ContainerInstanceArn' | awk -F/ '{print $2}' )",
      "sudo sed -i -e \"s/{container_instance_id}/$container_instance_id/g\" /etc/awslogs/awslogs.conf",
      "sudo service awslogs start",
      "sudo chkconfig awslogs on",
    ]

    connection {
      type = "ssh"
      user = "ec2-user"
      private_key = "${file("data-refinery-key.pem")}"
    }
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
