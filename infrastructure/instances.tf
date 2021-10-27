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
  template = file("workers-configuration/workers-instance-user-data.tpl.sh")

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

resource "aws_db_parameter_group" "postgres_parameters" {
  name = "postgres-parameters-${var.user}-${var.stage}"
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

  parameter {
    name = "log_checkpoints"
    value = false
  }

  parameter {
    name = "work_mem"
    value = 64000
  }

  # If you have a dedicated database server with 1GB or more of RAM, a
  # reasonable starting value for shared_buffers is 25% of the memory
  # in your system. There are some workloads where even larger
  # settings for shared_buffers are effective, but because PostgreSQL
  # also relies on the operating system cache, it is unlikely that an
  # allocation of more than 40% of RAM to shared_buffers will work
  # better than a smaller amount. Larger settings for shared_buffers
  # usually require a corresponding increase in max_wal_size, in order
  # to spread out the process of writing large quantities of new or
  # changed data over a longer period of time.
  parameter {
    name = "shared_buffers"
    # Note that the unit here is 8KB, so 8192 * .25 = 32768
    value = "{DBInstanceClassMemory/32768}"
  }

  # https://www.2ndquadrant.com/en/blog/basics-of-tuning-checkpoints/ says:
  # Let’s also discuss the other extreme – doing very frequent
  # checkpoints (say, every second or so). That would allow keeping
  # only tiny amount of WAL and the recovery would be very fast
  # (having to replay only the tiny WAL amount). But it would also
  # turn the asynchronous writes to data files into synchronous ones,
  # seriously impacting the users (e.g. increasing COMMIT latency,
  # reducing throughput).

  # Higher WAL values means a higher recovery time, but better general
  # performance. We want to optimize for general performance over
  # recovery time.
  parameter {
    name = "max_wal_size"
    # Note that the unit here is 1MB, so this is 4GB.
    value = 4096
  }

  parameter {
    name = "min_wal_size"
    # Note that the unit here is 1MB, so this is 1GB.
    value = 1024
  }

  parameter {
    name = "effective_cache_size"
    # Note that the unit here is 8KB, so 8192 * .5 = 16384
    value = "{DBInstanceClassMemory/16384}"
  }

  parameter {
    name = "wal_buffers"
    # Note that the unit here is 8KB so this is 16MB.
    value = 16384
  }

  parameter {
    name = "maintenance_work_mem"
    # The unit here is KB, so this is 1GB.
    value = 1048576
  }

  # From https://www.postgresql.org/docs/11/runtime-config-query.html
  # Storage that has a low random read cost relative to sequential,
  # e.g., solid-state drives, might also be better modeled with a
  # lower value for random_page_cost, e.g., 1.1.
  parameter {
    name = "random_page_cost"
    value = 1.1
  }

  # https://www.postgresql.org/docs/current/runtime-config-resource.html says
  # SSDs and other memory-based storage can often process many
  # concurrent requests, so the best value might be in the hundreds.
  parameter {
    name = "effective_io_concurrency"
    value = 200
  }

  parameter {
    name = "max_worker_processes"
    value = "{DBInstanceVCPU}"
  }

  parameter {
    name = "max_parallel_workers"
    value = "{DBInstanceVCPU}"
  }

  tags = var.default_tags
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

  apply_immediately = true

  db_subnet_group_name = aws_db_subnet_group.data_refinery.name
  parameter_group_name = aws_db_parameter_group.postgres_parameters.name

  # TF is broken, but we do want this protection in prod.
  # Related: https://github.com/hashicorp/terraform/issues/5417
  # Only the prod's bucket prefix is empty.
  skip_final_snapshot = var.stage == "prod" ? false : true
  final_snapshot_identifier = var.stage == "prod" ? "data-refinery-prod-snapshot" : "none"

  enabled_cloudwatch_logs_exports = ["postgresql", "upgrade"]

  vpc_security_group_ids = [aws_security_group.data_refinery_db.id]
  multi_az = true
  publicly_accessible = true

  backup_retention_period = var.stage == "prod" ? "7" : "0"

  tags = var.default_tags
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

  tags = merge(
    var.default_tags,
    {
      Name = "pg-bouncer-${var.user}-${var.stage}"
    }
  )

  root_block_device {
    volume_type = "gp2"
    volume_size = 100

    tags = var.default_tags
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
    instance_type = "m5.large.elasticsearch"
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

  tags = merge(
    var.default_tags,
    {
      Domain = "es-${var.user}-${var.stage}"
      Name = "es-${var.user}-${var.stage}"
    }
  )
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

  tags = merge(
    var.default_tags,
    {
      Name = "API Server 1 ${var.user}-${var.stage}"
    }
  )

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device {
    volume_type = "gp2"
    volume_size = 100

    tags = var.default_tags
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

  tags = merge(
    var.default_tags,
    {
      Name = "Foreman Server 1 ${var.user}-${var.stage}"
    }
  )

  # I think these are the defaults provided in terraform examples.
  # They should be removed or revisited.
  root_block_device {
    volume_type = "gp2"
    volume_size = 100

    tags = var.default_tags
  }
}

output "foreman_server_1_ip" {
  value = aws_instance.foreman_server_1.public_ip
}
