# The configuration contained in this file specifies the AWS instances
# we'll need. These include EC2 instances and an RDS instance.

data "aws_ami" "ubuntu" {
  most_recent = true

  filter {
    name   = "name"
    # This is a HVM, EBS backed SSD Ubuntu LTS AMI with Docker version 17.12.0 on it in the US,
    # the stock Ubuntu cloud image in the EU.
    values = ["${var.region == "us-east-1" ? "ubuntu-16.04-docker-17.12.0-ce-*" : "ubuntu/images/hvm-ssd/ubuntu-xenial-16.04-amd64-server-*" }"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

}

##
# Nomad Settings
##

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

##
# Nomad Cluster
##

# This instance will run the lead Nomad Server. We will want three
# Nomad Servers to prevent data loss, but one server needs to start
# first so that the others can be aware of its IP address to join it.
resource "aws_instance" "nomad_server_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "${var.nomad_server_instance_type}"
  availability_zone = "${var.region}b"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  disable_api_termination = "${var.stage == "prod" ? true : false}"

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
  instance_type = "${var.nomad_server_instance_type}"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  disable_api_termination = "${var.stage == "prod" ? true : false}"

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
}

output "nomad_server_2_ip" {
  value = "${aws_instance.nomad_server_2.public_ip}"
}

# This is another Nomad Server instance since it's recommended to have
# at least three to avoid data loss.
resource "aws_instance" "nomad_server_3" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "${var.nomad_server_instance_type}"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_worker.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_internet_gateway.data_refinery"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  disable_api_termination = "${var.stage == "prod" ? true : false}"

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
}

output "nomad_server_3_ip" {
  value = "${aws_instance.nomad_server_3.public_ip}"
}

# The Nomad Client needs to be aware of the Nomad Server's IP address,
# so we template it into its configuration.
data "template_file" "nomad_client_config" {
  template = "${file("nomad-configuration/client.tpl.hcl")}"

  vars {
    nomad_lead_server_ip = "${aws_instance.nomad_server_1.private_ip}"
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
    database_host = "${aws_db_instance.postgres_db.address}"
    database_port = "${var.database_hidden_port}"
    database_user = "${var.database_user}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
  }
}

data "template_file" "nomad_client_script_smasher_smusher" {
  template = "${file("nomad-configuration/client-smasher-instance-user-data.tpl.sh")}"

  vars {
    install_nomad_script = "${data.local_file.install_nomad_script.content}"
    nomad_client_config = "${data.template_file.nomad_client_config.rendered}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
    database_host = "${aws_db_instance.postgres_db.address}"
    database_port = "${var.database_hidden_port}"
    database_user = "${var.database_user}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
  }
}

##
# Autoscaling / Spot Fleets
##

resource "aws_spot_fleet_request" "cheap_ram" {
  iam_fleet_role      = "${aws_iam_role.data_refinery_spot_fleet.arn}"
  allocation_strategy = "diversified"
  target_capacity     = 100 # We're using RAM/10 here. (1000 is apparently too many for AWS.)
  valid_until         = "2021-11-04T20:44:20Z"
  fleet_type          = "maintain"

  # Instances won't be destroyed on Terraform destroy without this flag.
  # See https://github.com/hashicorp/terraform/issues/13859
  terminate_instances_with_expiration = true

  ##
  # Common / Depends On
  ##
  depends_on = [
            "aws_internet_gateway.data_refinery",
            "aws_instance.nomad_server_1",
            "aws_ebs_volume.data_refinery_ebs",
            "aws_instance.pg_bouncer"
  ]

  ##
  # x1.16xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "x1.16xlarge"
    weighted_capacity         = 10 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification x1.16xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }
  ##
  # x1.32xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "x1.32xlarge"
    weighted_capacity         = 20 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification x1.32xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }

  ##
  # x1e.8xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "x1e.8xlarge"
    weighted_capacity         = 10 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification x1e.8xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }

  ##
  # x1e.16xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "x1e.16xlarge"
    weighted_capacity         = 20 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification x1e.16xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }

  ##
  # x1e.32xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "x1e.32xlarge"
    weighted_capacity         = 40 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification x1e.32xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }

  ##
  # r5d.24xlarge
  ##
  launch_specification {

    # Client Specific
    instance_type             = "r5d.24xlarge"
    weighted_capacity         = 8 # Really, more like 1.4 # via https://aws.amazon.com/ec2/instance-types/
    spot_price                = "${var.spot_price}"
    ami                       = "${data.aws_ami.ubuntu.id}"
    iam_instance_profile_arn  = "${aws_iam_instance_profile.data_refinery_instance_profile.arn}"
    user_data                 = "${data.template_file.nomad_client_script_smusher.rendered}"
    vpc_security_group_ids    = ["${aws_security_group.data_refinery_worker.id}"]
    subnet_id                 = "${aws_subnet.data_refinery_1a.id}"
    availability_zone         = "${var.region}a"
    key_name = "${aws_key_pair.data_refinery.key_name}"

    root_block_device {
      volume_size = 20
      volume_type = "gp2"
    }

    tags {
        Name = "Spot Fleet Launch Specification r5d.24xlarge ${var.user}-${var.stage}"
        User = "${var.user}"
        Stage = "${var.stage}"
    }
  
  }

}

##
# Database
##

resource "aws_db_parameter_group" "postgres_parameters" {
  name = "postgres-parameters-${var.user}-${var.stage}"
  description = "Postgres Parameters ${var.user} ${var.stage}"
  family = "postgres9.6"

  parameter {
    name = "deadlock_timeout"
    value = "60000" # 60000ms = 60s
  }

  parameter {
    name = "statement_timeout"
    value = "60000" # 60000ms = 60s
  }
}

resource "aws_db_instance" "postgres_db" {
  identifier = "data-refinery-${var.user}-${var.stage}"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "9.6.6"
  instance_class = "db.${var.database_instance_type}"
  name = "data_refinery"
  port = "${var.database_hidden_port}"
  username = "${var.database_user}"
  password = "${var.database_password}"

  db_subnet_group_name = "${aws_db_subnet_group.data_refinery.name}"
  parameter_group_name = "${aws_db_parameter_group.postgres_parameters.name}"

  # TF is broken, but we do want this protection in prod.
  # Related: https://github.com/hashicorp/terraform/issues/5417
  # Only the prod's bucket prefix is empty.
  skip_final_snapshot = "${var.stage == "prod" ? false : true}"
  final_snapshot_identifier = "${var.stage == "prod" ? "data-refinery-prod-snapshot" : "none"}"

  vpc_security_group_ids = ["${aws_security_group.data_refinery_db.id}"]
  multi_az = true
  publicly_accessible = true

  backup_retention_period  = "${var.stage == "prod" ? "7" : "0"}"

}

resource "aws_instance" "pg_bouncer" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "${var.nomad_server_instance_type}"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_pg.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_db_instance.postgres_db"]
  key_name = "${aws_key_pair.data_refinery.key_name}"

  # Don't delete production servers by mistake!
  disable_api_termination = "${var.stage == "prod" ? true : false}"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = "${data.template_file.pg_bouncer_script_smusher.rendered}"

  tags = {
    Name = "pg-bouncer-${var.user}-${var.stage}"
  }

  # Nomad server requirements can be found here:
  # https://www.nomadproject.io/guides/cluster/requirements.html
  # However I do not think that these accurately reflect those requirements.
  # I think these are the defaults provided in terraform examples.
  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }
}

data "template_file" "pg_bouncer_script_smusher" {
  template = "${file("nomad-configuration/pg-bouncer-instance-user-data.tpl.sh")}"

  vars {
    database_host = "${aws_db_instance.postgres_db.address}"
    database_user = "${var.database_user}"
    database_port = "${var.database_hidden_port}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
    listen_port = "${var.database_port}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
  }
}

##
# API Webserver
##

data "local_file" "api_nginx_config" {
  filename = "api-configuration/nginx_config.conf"
}

data "local_file" "api_environment" {
  filename = "api-configuration/environment"
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad_lead_server_script_smusher} but for the Nginx/API.
data "template_file" "api_server_script_smusher" {
  template = "${file("api-configuration/api-server-instance-user-data.tpl.sh")}"

  vars {
    nginx_config = "${data.local_file.api_nginx_config.content}"
    api_environment = "${data.local_file.api_environment.content}"
    dockerhub_repo = "${var.dockerhub_repo}"
    api_docker_image = "${var.api_docker_image}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
    database_host = "${aws_instance.pg_bouncer.private_ip}"
    database_user = "${var.database_user}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
    log_group = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
    log_stream = "${aws_cloudwatch_log_stream.log_stream_api.name}"
  }
}

resource "aws_instance" "api_server_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "${var.api_instance_type}"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_api.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = [
    "aws_db_instance.postgres_db",
    "aws_instance.pg_bouncer",
    "aws_security_group_rule.data_refinery_api_http",
    "aws_security_group_rule.data_refinery_api_outbound"
  ]
  user_data = "${data.template_file.api_server_script_smusher.rendered}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "API Server 1 ${var.user}-${var.stage}"
  }

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }
}

output "api_server_1_ip" {
  value = "${aws_instance.api_server_1.public_ip}"
}

##
# Foreman Server
##

data "local_file" "foreman_environment" {
  filename = "foreman-configuration/environment"
}

# This script smusher serves a similar purpose to
# ${data.template_file.nomad_lead_server_script_smusher} but for the Foreman.
data "template_file" "foreman_server_script_smusher" {
  template = "${file("foreman-configuration/foreman-server-instance-user-data.tpl.sh")}"

  vars {
    foreman_environment = "${data.local_file.foreman_environment.content}"
    dockerhub_repo = "${var.dockerhub_repo}"
    foreman_docker_image = "${var.foreman_docker_image}"
    user = "${var.user}"
    stage = "${var.stage}"
    region = "${var.region}"
    nomad_lead_server_ip = "${aws_instance.nomad_server_1.private_ip}"
    database_host = "${aws_instance.pg_bouncer.private_ip}"
    database_user = "${var.database_user}"
    database_password = "${var.database_password}"
    database_name = "${aws_db_instance.postgres_db.name}"
    log_group = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
  }
}

resource "aws_instance" "foreman_server_1" {
  ami = "${data.aws_ami.ubuntu.id}"
  instance_type = "${var.foreman_instance_type}"
  availability_zone = "${var.region}a"
  vpc_security_group_ids = ["${aws_security_group.data_refinery_foreman.id}"]
  iam_instance_profile = "${aws_iam_instance_profile.data_refinery_instance_profile.name}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  depends_on = ["aws_db_instance.postgres_db", "aws_instance.pg_bouncer"]
  user_data = "${data.template_file.foreman_server_script_smusher.rendered}"
  key_name = "${aws_key_pair.data_refinery.key_name}"

  tags = {
    Name = "Foreman Server 1 ${var.user}-${var.stage}"
  }

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device = {
    volume_type = "gp2"
    volume_size = 100
  }
}

output "foreman_server_1_ip" {
  value = "${aws_instance.foreman_server_1.public_ip}"
}
