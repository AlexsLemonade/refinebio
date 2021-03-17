# The configuration contained in this file specifies the AWS instances
# we'll need. These include EC2 instances and an RDS instance.

data "aws_ami" "ubuntu" {
  most_recent = true

  owners = ["589864003899"]

  filter {
    name = "name"
    values = ["ccdl-ubuntu-18.04-*"]
  }
}

# This script smusher exists in order to be able to circumvent a
# limitation of AWS which is that you get one script and one script
# only to set up the instance when it boots up. Because there is only
# one script you cannot place additional files your script may need
# onto the instance. Therefore this script smusher templates the files
# the instance-user-data.sh script needs into it, so that once it
# makes its way onto the instance it can spit them back out onto the
# disk.
data "template_file" "worker_script_smusher" {
  template = file("workers-configuration/client-instance-user-data.tpl.sh")

  vars = {
    user = var.user
    stage = var.stage
    region = var.region
    database_host = aws_instance.pg_bouncer.private_ip
    database_port = var.database_port
    database_user = var.database_user
    database_password = var.database_password
    database_name = aws_db_instance.postgres_db.name
  }
}


##
# Database
##

# Temporary. This will go away once we can delete things.
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

resource "aws_db_parameter_group" "postgres_parameter_group" {
  name = "postgres-parameter-group-${var.user}-${var.stage}"
  description = "Postgres Parameters ${var.user} ${var.stage}"
  family = "postgres11"

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
  engine_version = "11.1"
  allow_major_version_upgrade = true
  auto_minor_version_upgrade = false
  instance_class = "db.${var.database_instance_type}"
  name = "data_refinery"
  port = var.database_hidden_port
  username = var.database_user
  password = var.database_password

  db_subnet_group_name = aws_db_subnet_group.data_refinery.name
  parameter_group_name = aws_db_parameter_group.postgres_parameter_group.name

  # TF is broken, but we do want this protection in prod.
  # Related: https://github.com/hashicorp/terraform/issues/5417
  # Only the prod's bucket prefix is empty.
  skip_final_snapshot = var.stage == "prod" ? false : true
  final_snapshot_identifier = var.stage == "prod" ? "data-refinery-prod-snapshot" : "none"

  vpc_security_group_ids = [aws_security_group.data_refinery_db.id]
  multi_az = true
  publicly_accessible = true

  backup_retention_period = var.stage == "prod" ? "7" : "0"
}

resource "aws_instance" "pg_bouncer" {
  ami = data.aws_ami.ubuntu.id
  instance_type = var.pg_bouncer_instance_type
  availability_zone = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_pg.id]
  iam_instance_profile = aws_iam_instance_profile.data_refinery_api.name
  subnet_id = aws_subnet.data_refinery_1a.id
  depends_on = [aws_db_instance.postgres_db]
  key_name = aws_key_pair.data_refinery.key_name
  disable_api_termination = "false"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = data.template_file.pg_bouncer_script_smusher.rendered

  tags = {
    Name = "pg-bouncer-${var.user}-${var.stage}"
  }

  root_block_device {
    volume_type = "gp2"
    volume_size = 100
  }
}

data "template_file" "pg_bouncer_script_smusher" {
  template = file("workers-configuration/pg-bouncer-instance-user-data.tpl.sh")

  vars = {
    database_host = aws_db_instance.postgres_db.address
    database_user = var.database_user
    database_port = var.database_hidden_port
    database_password = var.database_password
    database_name = aws_db_instance.postgres_db.name
    listen_port = var.database_port
    user = var.user
    stage = var.stage
    region = var.region
  }
}

##
# ElasticSearch
##

# This only needs to be run one time per account.
# Run it zero times, it won't work. Run it more than one time, it won't work.
# TF needs away to manage these.

# Related: https://github.com/terraform-providers/terraform-provider-aws/issues/5218
# Related: https://github.com/cloudposse/terraform-aws-elasticsearch/issues/5
# resource "aws_iam_service_linked_role" "es" {
#   aws_service_name = "es.amazonaws.com"
# }

data "aws_caller_identity" "current" {
}

resource "aws_elasticsearch_domain" "es" {
  domain_name = "es-${var.user}-${var.stage}"
  elasticsearch_version = "6.3"

  advanced_options = {
    "indices.query.bool.max_clause_count" = 16384
  }

  # TODO: Figure out the power/cost balance of this type.
  # Prices are here: https://aws.amazon.com/elasticsearch-service/pricing/
  cluster_config {
    instance_type = "t2.medium.elasticsearch"
  }

  vpc_options {
    subnet_ids = [
      aws_subnet.data_refinery_1a.id,
    ]
    security_group_ids = [
      aws_security_group.data_refinery_es.id,
    ]
  }

  ebs_options {
    ebs_enabled = true

    # This depends on the instance type, else you'll get this error:
    #   * aws_elasticsearch_domain.es: LimitExceededException: Volume size must be between 10 and 35 for t2.medium.elasticsearch instance type and elasticsearch version 6.3
    volume_size = 10
  }

  access_policies = <<CONFIG
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "AWS": "*"
      },
      "Action": "es:*",
      "Resource": "arn:aws:es:${var.region}:${data.aws_caller_identity.current.account_id}:domain/es-${var.user}-${var.stage}/*"
    }
  ]
}

CONFIG


  snapshot_options {
    automated_snapshot_start_hour = 23
  }

  tags = {
    Domain = "es-${var.user}-${var.stage}"
    Name = "es-${var.user}-${var.stage}"
  }
}

output "elasticsearch_endpoint" {
  value = aws_elasticsearch_domain.es.endpoint
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
# ${data.template_file.worker_script_smusher} but for the Nginx/API.
data "template_file" "api_server_script_smusher" {
  template = file("api-configuration/api-server-instance-user-data.tpl.sh")

  vars = {
    nginx_config = data.local_file.api_nginx_config.content
    api_environment = data.local_file.api_environment.content
    dockerhub_repo = var.dockerhub_repo
    api_docker_image = var.api_docker_image
    user = var.user
    stage = var.stage
    region = var.region
    database_host = aws_instance.pg_bouncer.private_ip
    database_user = var.database_user
    database_password = var.database_password
    database_name = aws_db_instance.postgres_db.name
    elasticsearch_host = aws_elasticsearch_domain.es.endpoint
    elasticsearch_port = "80" # AWS doesn't support the data transfer protocol on 9200 >:[
    log_group = aws_cloudwatch_log_group.data_refinery_log_group.name
    log_stream = aws_cloudwatch_log_stream.log_stream_api.name
  }

  depends_on = [
    aws_db_instance.postgres_db,
    aws_elasticsearch_domain.es,
    aws_instance.pg_bouncer,
    aws_security_group_rule.data_refinery_api_http,
    aws_security_group_rule.data_refinery_api_outbound,
  ]
}

resource "aws_instance" "api_server_1" {
  ami = data.aws_ami.ubuntu.id
  instance_type = var.api_instance_type
  availability_zone = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_api.id]
  iam_instance_profile = aws_iam_instance_profile.data_refinery_api.name
  subnet_id = aws_subnet.data_refinery_1a.id
  depends_on = [
    aws_db_instance.postgres_db,
    aws_elasticsearch_domain.es,
    aws_instance.pg_bouncer,
    aws_security_group_rule.data_refinery_api_http,
    aws_security_group_rule.data_refinery_api_outbound,
  ]
  user_data = data.template_file.api_server_script_smusher.rendered
  key_name = aws_key_pair.data_refinery.key_name

  tags = {
    Name = "API Server 1 ${var.user}-${var.stage}"
  }

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device {
    volume_type = "gp2"
    volume_size = 100
  }
}

output "api_server_1_ip" {
  value = aws_instance.api_server_1.public_ip
}

##
# Foreman Server
##

data "local_file" "foreman_environment" {
  filename = "foreman-configuration/environment"
}

# This script smusher serves a similar purpose to
# ${data.template_file.worker_script_smusher} but for the Foreman.
data "template_file" "foreman_server_script_smusher" {
  template = file(
    "foreman-configuration/foreman-server-instance-user-data.tpl.sh",
  )

  vars = {
    foreman_environment = data.local_file.foreman_environment.content
    dockerhub_repo = var.dockerhub_repo
    foreman_docker_image = var.foreman_docker_image
    user = var.user
    stage = var.stage
    region = var.region
    database_host = aws_instance.pg_bouncer.private_ip
    database_user = var.database_user
    database_password = var.database_password
    database_name = aws_db_instance.postgres_db.name
    elasticsearch_host = aws_elasticsearch_domain.es.endpoint
    elasticsearch_port = var.elasticsearch_port
    log_group = aws_cloudwatch_log_group.data_refinery_log_group.name
    aws_access_key_id = aws_iam_access_key.data_refinery_user_client_key.id
    aws_secret_access_key = aws_iam_access_key.data_refinery_user_client_key.secret
  }
}

resource "aws_instance" "foreman_server_1" {
  ami = data.aws_ami.ubuntu.id
  instance_type = var.foreman_instance_type
  availability_zone = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_foreman.id]
  iam_instance_profile = aws_iam_instance_profile.data_refinery_foreman.name
  subnet_id = aws_subnet.data_refinery_1a.id
  depends_on = [
    aws_db_instance.postgres_db,
    aws_instance.pg_bouncer,
    aws_elasticsearch_domain.es,
  ]
  user_data = data.template_file.foreman_server_script_smusher.rendered
  key_name = aws_key_pair.data_refinery.key_name

  tags = {
    Name = "Foreman Server 1 ${var.user}-${var.stage}"
  }

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device {
    volume_type = "gp2"
    volume_size = 100
  }
}

output "foreman_server_1_ip" {
  value = aws_instance.foreman_server_1.public_ip
}
