provider "aws" {
  version = "2.70.0"
  region = "us-east-1"
}

data "aws_vpc" "ccdl_dev_vpc" {
  filter {
    name = "tag:Name"
    values = ["ccdl-dev-vpc"]
  }
}

data "aws_subnet" "ccdl_dev_subnet" {
  filter {
    name = "tag:Name"
    values = ["ccdl-dev-vpc-subnet"]
  }
}

resource "aws_instance" "ubuntu-ami-template-instance" {
  # this is the base ubuntu AMI on AWS as of 2020-07-10
  ami = "ami-0ac80df6eff0e70b5"

  instance_type = "t2.micro"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = data.template_file.ubuntu_instance_user_data.rendered

  subnet_id = data.aws_subnet.ccdl_dev_subnet.id
  associate_public_ip_address = true
  key_name = "data-refinery-key-circleci-prod"
  vpc_security_group_ids = [aws_security_group.ami_template_instance.id]

  tags = {
    Name = "Ubuntu AMI Template Instance"
  }
}

resource "aws_instance" "ecs-ami-template-instance" {
  # This is the Amazon ECS-optimized Amazon Linux 2 AMI as of 2020-12-17
  ami = "ami-0f06fc190dd71269e"

  instance_type = "t2.micro"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = data.template_file.ecs_instance_user_data.rendered

  subnet_id = data.aws_subnet.ccdl_dev_subnet.id
  associate_public_ip_address = true
  key_name = "data-refinery-key-circleci-prod"
  vpc_security_group_ids = [aws_security_group.ami_template_instance.id]

  tags = {
    Name = "ECS AMI Template Instance"
  }
}
