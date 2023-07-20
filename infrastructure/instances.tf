# The configuration contained in this file specifies the AWS instances
# we'll need. These include EC2 instances and an RDS instance.

data "aws_ami" "ubuntu" {
  most_recent = true

  owners = ["589864003899"]

  filter {
    name   = "name"
    values = ["ccdl-ubuntu-18.04-*"]
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
  domain_name           = "es-${var.user}-${var.stage}"
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
      Name   = "es-${var.user}-${var.stage}"
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

resource "aws_instance" "api_server_1" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.api_instance_type
  availability_zone      = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_api.id]
  iam_instance_profile   = aws_iam_instance_profile.data_refinery_api.name
  subnet_id              = aws_subnet.data_refinery_1a.id
  depends_on = [
    aws_db_instance.postgres_db,
    aws_elasticsearch_domain.es,
    aws_instance.pg_bouncer,
    aws_s3_bucket.data_refinery_cert_bucket,
    aws_security_group_rule.data_refinery_api_http,
    aws_security_group_rule.data_refinery_api_outbound,
  ]

  user_data = templatefile("api-configuration/api-server-instance-user-data.tpl.sh",
    {
      api_docker_image          = var.api_docker_image
      api_environment           = data.local_file.api_environment.content
      data_refinery_cert_bucket = aws_s3_bucket.data_refinery_cert_bucket.id
      database_host             = aws_instance.pg_bouncer.private_ip
      database_name             = aws_db_instance.postgres_db.db_name
      database_password         = var.database_password
      database_user             = var.database_user
      dockerhub_repo            = var.dockerhub_repo
      elasticsearch_host        = aws_elasticsearch_domain.es.endpoint
      elasticsearch_port        = "80" # AWS doesn't support the data transfer protocol on 9200 >:[
      log_group                 = aws_cloudwatch_log_group.data_refinery_log_group.name
      log_stream                = aws_cloudwatch_log_stream.log_stream_api.name
      nginx_config              = data.local_file.api_nginx_config.content
      region                    = var.region
      stage                     = var.stage
      user                      = var.user
    }
  )
  key_name  = aws_key_pair.data_refinery.key_name

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

resource "aws_instance" "foreman_server_1" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.foreman_instance_type
  availability_zone      = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_foreman.id]
  iam_instance_profile   = aws_iam_instance_profile.data_refinery_foreman.name
  subnet_id              = aws_subnet.data_refinery_1a.id

  depends_on = [
    aws_db_instance.postgres_db,
    aws_elasticsearch_domain.es,
    aws_instance.pg_bouncer,
    aws_s3_bucket.data_refinery_cert_bucket,
    aws_security_group_rule.data_refinery_api_http,
    aws_security_group_rule.data_refinery_api_outbound,
  ]

  user_data = templatefile("foreman-configuration/foreman-server-instance-user-data.tpl.sh",
    {
      accession_gathering_job_run_day = var.accession_gathering_job_run_day
      database_host = aws_instance.pg_bouncer.private_ip
      database_name = aws_db_instance.postgres_db.db_name
      database_password = var.database_password
      database_user = var.database_user
      dockerhub_repo = var.dockerhub_repo
      elasticsearch_host = aws_elasticsearch_domain.es.endpoint
      elasticsearch_port = var.elasticsearch_port
      foreman_docker_image = var.foreman_docker_image
      foreman_environment = data.local_file.foreman_environment.content
      log_group = aws_cloudwatch_log_group.data_refinery_log_group.name
      region = var.region
      stage = var.stage
      user = var.user
    }
  )
  key_name  = aws_key_pair.data_refinery.key_name

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
