# The configuration contained in this file specifies the AWS instances
# we'll need. These include EC2 instances and an RDS instance.

data "aws_ami" "ubuntu" {
  most_recent = true

  filter {
    name   = "name"
    # This is a HVM, EBS backed SSD Ubuntu LTS AMI with Docker version 17.12.0 on it.
    values = ["ubuntu-16.04-docker-17.12.0-ce-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }
}

# We need to transfer the Nomad Job Specifications onto the instance
# so we can register them. However we use these job specifications to
# store environment variables, so they are on Github as templates. To
# format them with the production environment variables, we need to
# run the following command. By doing so with a `null_resource` we're
# able to avoid having a prerequisite step to running `terraform
# apply`.
resource "null_resource" "format-nomad-job-specs" {
  provisioner "local-exec" {
    command = "cd .. && REGION=${var.region} USER=${var.user} STAGE=${var.stage} ./workers/format_nomad_with_env.sh -e prod -o $(pwd)/infrastructure/nomad-job-specs/"
  }
}

# This script installs Nomad.
data "local_file" "install-nomad-script" {
  filename = "../install_nomad.sh"
}

# This is the configuration for the Nomad Server.
data "local_file" "nomad-lead-server-config" {
  filename = "nomad-configuration/lead_server.hcl"
}

# This is a Nomad Job Specification file built by ${null_resource.format-nomad-job-specs}.
data "local_file" "downloader-job-spec" {
  filename = "nomad-job-specs/downloader.nomad"
}

# This is another Nomad Job Specification file built by ${null_resource.format-nomad-job-specs}.
data "local_file" "processor-job-spec" {
  filename = "nomad-job-specs/processor.nomad"
}

# This script smusher exists in order to be able to circumvent a
# limitation of AWS which is that you get one script and one script
# only to set up the instance when it boots up. Because there is only
# one script you cannot place additional files your script may need
# onto the instance. Therefore this script smusher templates the files
# the instance-user-data.sh script needs into it, so that once it
# makes its way onto the instance it can spit them back out onto the
# disk.
data "template_file" "nomad-lead-server-script-smusher" {
  template = "${file("nomad-configuration/lead-server-instance-user-data.tpl.sh")}"

  vars {
    downloader_job_spec = "${data.local_file.downloader-job-spec.content}"
    processor_job_spec = "${data.local_file.processor-job-spec.content}"
    install_nomad_script = "${data.local_file.install-nomad-script.content}"
    nomad_server_config = "${data.local_file.nomad-lead-server-config.content}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
  }
}

# This instance will run the lead Nomad Server. We will want three
# Nomad Servers to prevent data loss, but one server needs to start
# first so that the others can be aware of its IP address to join it.
resource "aws_instance" "nomad-server-1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad-lead-server-script-smusher.rendered}"

  tags = {
    Name = "nomad-server-1-${var.user}-${var.stage}"
  }

  # Nomad server requirements can be found here:
  # https://www.nomadproject.io/guides/cluster/requirements.html
  # However I do not think that these accurately reflect those requirements.
  # I think these are the defaults provided in terraform examples.
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

output "nomad-server-1-ip" {
  value = "${aws_instance.nomad-server-1.public_ip}"
}

# The other Nomad Servers needs to be aware of the Lead Nomad Server's
# IP address, so we template it into their configurations.
data "template_file" "nomad-server-config" {
  template = "${file("nomad-configuration/server.tpl.hcl")}"

  vars {
    nomad_lead_server_ip = "${aws_instance.nomad-server-1.private_ip}"
  }
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad-lead-server-script-smusher} but for the other
# Nomad Servers. See that resource's definition for more information.
data "template_file" "nomad-server-script-smusher" {
  template = "${file("nomad-configuration/server-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install-nomad-script.content}"
    nomad_server_config = "${data.template_file.nomad-server-config.rendered}"
  }
}

# This is another Nomad Server instance since it's recommended to have
# at least three to avoid data loss.
resource "aws_instance" "nomad-server-2" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "us-east-1a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad-server-script-smusher.rendered}"

  tags = {
    Name = "nomad-server-2-${var.user}-${var.stage}"
  }

  # Nomad server requirements can be found here:
  # https://www.nomadproject.io/guides/cluster/requirements.html
  # However I do not think that these accurately reflect those requirements.
  # I think these are the defaults provided in terraform examples.
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

output "nomad-server-2-ip" {
  value = "${aws_instance.nomad-server-2.public_ip}"
}

# This is another Nomad Server instance since it's recommended to have
# at least three to avoid data loss.
resource "aws_instance" "nomad-server-3" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "us-east-1a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad-server-script-smusher.rendered}"

  tags = {
    Name = "nomad-server-3-${var.user}-${var.stage}"
  }

  # Nomad server requirements can be found here:
  # https://www.nomadproject.io/guides/cluster/requirements.html
  # However I do not think that these accurately reflect those requirements.
  # I think these are the defaults provided in terraform examples.
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

output "nomad-server-3-ip" {
  value = "${aws_instance.nomad-server-3.public_ip}"
}

# The Nomad Client needs to be aware of the Nomad Server's IP address,
# so we template it into its configuration.
data "template_file" "nomad-client-config" {
  template = "${file("nomad-configuration/client.tpl.hcl")}"

  vars {
    nomad_server_address = "${aws_instance.nomad-server-1.private_ip}"
  }
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad-lead-server-script-smusher} but for the Nomad
# Client. See that resource's definition for more information.
data "template_file" "nomad-client-script-smusher" {
  template = "${file("nomad-configuration/client-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install-nomad-script.content}"
    nomad_client_config = "${data.template_file.nomad-client-config.rendered}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
  }
}

# This instance will run the Nomad Client. We will eventually want
# these to be spot instances, but it is much simpler operationally to
# do it using a normal ec2 instance.
resource "aws_instance" "nomad-client-1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.xlarge"
  availability_zone = "${var.region}b"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.ecs_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  depends_on = ["aws_internet_gateway.data_refinery", "aws_instance.nomad-server-1"]
  user_data = "${data.template_file.nomad-client-script-smusher.rendered}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "nomad-client-1-${var.user}-${var.stage}"
  }

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
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



output "nomad-client-ip" {
  value = "${aws_instance.nomad-client-1.public_ip}"
}

variable "database_password" {}

resource "aws_db_instance" "postgres-db" {
  identifier = "data-refinery-${var.user}-${var.stage}"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.5.4"
  instance_class = "db.t2.micro"
  name = "data_refinery_${var.user}_${var.stage}"
  username = "data_refinery_user"
  password = "${var.database_password}"
  db_subnet_group_name = "${aws_db_subnet_group.data_refinery.name}"
  skip_final_snapshot = true
  vpc_security_group_ids = ["${aws_security_group.data_refinery_db.id}"]
  multi_az = true
  publicly_accessible = true
}
