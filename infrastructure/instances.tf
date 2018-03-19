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

# This script installs Nomad.
data "local_file" "install_nomad_script" {
  filename = "../install_nomad.sh"
}

# This is the configuration for the Nomad Server.
data "local_file" "nomad_lead_server_config" {
  filename = "nomad-configuration/lead_server.hcl"
}

# This script smusher exists in order to be able to circumvent a
# limitation of AWS which is that you get one script and one script
# only to set up the instance when it boots up. Because there is only
# one script you cannot place additional files your script may need
# onto the instance. Therefore this script smusher templates the files
# the instance-user-data.sh script needs into it, so that once it
# makes its way onto the instance it can spit them back out onto the
# disk.
data "template_file" "nomad_lead_server_script_smusher" {
  template = "${file("nomad-configuration/lead-server-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install_nomad_script.content}"
    nomad_server_config = "${data.local_file.nomad_lead_server_config.content}"
    server_number = 1
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
  }
}

# This instance will run the lead Nomad Server. We will want three
# Nomad Servers to prevent data loss, but one server needs to start
# first so that the others can be aware of its IP address to join it.
resource "aws_instance" "nomad_server_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad_lead_server_script_smusher.rendered}"

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

output "nomad_server_1_ip" {
  value = "${aws_instance.nomad_server_1.public_ip}"
}

# The other Nomad Servers needs to be aware of the Lead Nomad Server's
# IP address, so we template it into their configurations.
data "template_file" "nomad_server_config" {
  template = "${file("nomad-configuration/server.tpl.hcl")}"

  vars {
    nomad_lead_server_ip = "${aws_instance.nomad_server_1.private_ip}"
  }
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad_lead_server_script_smusher} but for the other
# Nomad Servers. See that resource's definition for more information.
data "template_file" "nomad_server_script_smusher" {
  template = "${file("nomad-configuration/server-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install_nomad_script.content}"
    nomad_server_config = "${data.template_file.nomad_server_config.rendered}"
    server_number = 1
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
  }
}

# This is another Nomad Server instance since it's recommended to have
# at least three to avoid data loss.
resource "aws_instance" "nomad_server_2" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "${var.region}b"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad_server_script_smusher.rendered}"

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

output "nomad_server_2_ip" {
  value = "${aws_instance.nomad_server_2.public_ip}"
}

# This is another Nomad Server instance since it's recommended to have
# at least three to avoid data loss.
resource "aws_instance" "nomad_server_3" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.small"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.nomad_server_script_smusher.rendered}"

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

output "nomad_server_3_ip" {
  value = "${aws_instance.nomad_server_3.public_ip}"
}

# The Nomad Client needs to be aware of the Nomad Server's IP address,
# so we template it into its configuration.
data "template_file" "nomad_client_config" {
  template = "${file("nomad-configuration/client.tpl.hcl")}"

  vars {
    nomad_server_address = "${aws_instance.nomad_server_1.private_ip}"
  }
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad_lead_server_script_smusher} but for the Nomad
# Client. See that resource's definition for more information.
data "template_file" "nomad_client_script_smusher" {
  template = "${file("nomad-configuration/client-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install_nomad_script.content}"
    nomad_client_config = "${data.template_file.nomad_client_config.rendered}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
    file_system_id = "${aws_efs_file_system.data_refinery_efs.id}"
    database_host = "${aws_db_instance.postgres_db.address}"
    database_user = "${var.database_user}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
  }
}

# This instance will run the Nomad Client. We will eventually want
# these to be spot instances, but it is much simpler operationally to
# do it using a normal ec2 instance.
resource "aws_instance" "nomad_client_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "t2.xlarge"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery", "aws_instance.nomad_server_1"]
  user_data = "${data.template_file.nomad_client_script_smusher.rendered}"
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

output "nomad_client_ip" {
  value = "${aws_instance.nomad_client_1.public_ip}"
}

resource "aws_db_instance" "postgres_db" {
  identifier = "data-refinery-${var.user}-${var.stage}"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.6.6"
  instance_class = "db.t2.micro"
  name = "data_refinery"
  username = "${var.database_user}"
  password = "${var.database_password}"
  db_subnet_group_name = "${aws_db_subnet_group.data_refinery.name}"
  # We probably actually want to keep this, but TF is broken here.
  # Related: https://github.com/hashicorp/terraform/issues/5417
  skip_final_snapshot = true
  vpc_security_group_ids = ["${aws_security_group.data_refinery_db.id}"]
  multi_az = true
  publicly_accessible = true
}
