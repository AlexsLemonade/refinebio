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

variable "database_password" {}

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
