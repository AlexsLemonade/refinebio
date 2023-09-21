# Variables
# Only lowercase alphanumeric characters and hyphens allowed!
# In CI, this can be set with TF_VAR_user, TF_VAR_stage, etc.

variable "default_tags" {
  default = {
    team    = "engineering"
    project = "refine.bio"
  }
  description = "Default resource tags"
  type        = map(string)
}

# Annoyingly, TF can't have computed variables (${env.USER})
# as default values for variables. So `TF_VAR_user=rjones tf plan`, etc.
# Also, don't use any non-alphanumeric characters here or RDS will whinge.
variable "user" {
  default = "myusername"
}

variable "stage" {
  default = "dev"
}

variable "environment" {
  default = "dev"
}

variable "static_bucket_prefix" {
  default = "staging"
}

variable "static_bucket_root" {
  default = ".refine.bio"
}

variable "region" {
  default = "us-east-1"
}

variable "host_ip" {
  # This will be overwritten.
  default = "127.0.0.1"
}

variable "ssh_public_key" {
  default = "MISSING_VALUE"
}

variable "database_user" {
  default = "drpostgresuser"
}

variable "database_password" {
  # This will be overwritten by the password in terraform.tfvars.
  # It's kept there so it's secret.
  default = "drpostgrespassword"
}

variable "django_secret_key" {
  # TODO: generate a new one of these and store it in terraform.tfvars as well.
  default = "NtG1bxZU115GThwrLuAJe0PhTVN9hJ4P"
}

variable "django_debug" {
  default = "False"
}

variable "database_port" {
  default = "5432"
}

# This is the port the RDS instance uses, but it is hidden from most
# of the application by PGBouncer.
variable "database_hidden_port" {
  default = "5430"
}

variable "database_timeout" {
  default = "30"
}

variable "database_instance_type" {
  default = "t2.micro"
}

variable "running_in_cloud" {
  default = "True"
}

variable "log_level" {
  default = "WARN"
}

variable "dockerhub_repo" {
  default = "ccdlstaging"
}

variable "downloaders_docker_image" {
  default = "dr_downloaders:latest"
}

variable "transcriptome_docker_image" {
  default = "dr_transcriptome:latest"
}

variable "salmon_docker_image" {
  default = "dr_salmon:latest"
}

variable "smasher_docker_image" {
  default = "dr_smasher:latest"
}

variable "affymetrix_docker_image" {
  default = "dr_affymetrix:latest"
}

variable "illumina_docker_image" {
  default = "dr_illumina:latest"
}

variable "no_op_docker_image" {
  default = "dr_no_op:latest"
}

variable "foreman_docker_image" {
  default = "dr_foreman:latest"
}

variable "use_s3" {
  default = "True"
}

variable "local_root_dir" {
  default = "/home/user/data_store"
}

# Instance types / ASG
variable "pg_bouncer_instance_type" {
  default = "t2.medium"
}

variable "smasher_instance_type" {
  # 128GiB Memory, smasher and compendia jobs need 30.
  # RNA-seq compendia needs 131gb
  # Enough for one smasher job at a time:
  default = "m5.2xlarge"
  # Appropriate for most compendia.
  # default = "m5.8xlarge"
  # Required for human and mouse quantpendia.
  # default = "m5.16xlarge"

  # 976GiB Memory, smasher and compendia jobs need 900.
  # Required for human and mouse compendia
  # default = "x1.16xlarge"
}

variable "spot_fleet_capacity" {
  default = "0"
}

# AWS Secrets manager won't allow these to be empty, so we make them
# "None" and check for that string in the settings,
variable "sentry_dsn" {
  default = "None"
}

# API
variable "api_docker_image" {
  default = "dr_api:latest"
}

variable "api_instance_type" {
  default = "t2.large"
}

variable "foreman_instance_type" {
  default = "m5.2xlarge"
}

variable "worker_ami" {
  default = "ami-0f9b369e4ad339b15"
}

variable "smasher_volume_size_in_gb" {
  # 500 is the smallest for ST1s.
  default = "500"
}

variable "num_workers" {
  # For now err on the side of inexpensive
  default = 1
}

variable "batch_use_on_demand_instances" {
  type    = bool
  default = false
}

variable "max_jobs_per_node" {
  default = 250
}

variable "max_downloader_jobs_per_node" {
  default = 200
}

variable "elasticsearch_port" {
  default = "80"
}

variable "slack_webhook_url" {
  # Only necessary for TF, but will be overwritten.
  default = "DEFAULT"
}

variable "full_stack" {
  default = false
}

variable "processing_compendia" {
  default = true
}

variable "accession_gathering_job_run_day" {
  default = "SAT"
}

variable "max_accessions_gathered_per_run" {
  default = 0
}

# Output our production environment variables.
output "environment_variables" {
  value = [
    {
      name  = "AWS_REGION"
      value = var.region
    },
    {
      name  = "USER"
      value = var.user
    },
    {
      name  = "STAGE"
      value = var.stage
    },
    {
      name  = "SSH_PUBLIC_KEY"
      value = var.ssh_public_key
    },
    {
      name  = "WORKER_ROLE_ARN"
      value = aws_iam_role.data_refinery_worker.arn
    },
    {
      name  = "DJANGO_DEBUG"
      value = var.django_debug
    },
    {
      name  = "DJANGO_SECRET_KEY"
      value = var.django_secret_key
    },
    {
      name  = "DJANGO_SECRET_KEY_ARN"
      value = aws_secretsmanager_secret.django_secret_key.arn
    },
    {
      name  = "DATABASE_NAME"
      value = aws_db_instance.postgres_db.db_name
    },
    {
      name  = "DATABASE_HOST"
      value = aws_instance.pg_bouncer.private_ip
    },
    {
      # We want to use the private IP from everywhere except wherever
      # deployment is happening, so we also need to expose this IP.
      name  = "DATABASE_PUBLIC_HOST"
      value = aws_instance.pg_bouncer.public_ip
    },
    {
      name  = "RDS_HOST"
      value = aws_db_instance.postgres_db.address
    },
    {
      name  = "DATABASE_USER"
      value = var.database_user
    },
    {
      name  = "DATABASE_PASSWORD"
      value = var.database_password
    },
    {
      name  = "DATABASE_PASSWORD_ARN"
      value = aws_secretsmanager_secret.database_password.arn
    },
    {
      name  = "DATABASE_PORT"
      value = var.database_port
    },
    {
      name  = "DATABASE_HIDDEN_PORT"
      value = var.database_hidden_port
    },
    {
      name  = "DATABASE_TIMEOUT"
      value = var.database_timeout
    },
    {
      name  = "API_HOST"
      value = aws_eip.data_refinery_api_ip.public_ip
    },
    {
      name  = "ELASTICSEARCH_HOST"
      value = aws_elasticsearch_domain.es.endpoint
    },
    {
      name  = "ELASTICSEARCH_PORT"
      value = var.elasticsearch_port
    },
    {
      name  = "RUNNING_IN_CLOUD"
      value = var.running_in_cloud
    },
    {
      name  = "LOG_LEVEL"
      value = var.log_level
    },
    {
      name  = "USE_S3"
      value = var.use_s3
    },
    {
      name  = "SENTRY_DSN"
      value = var.sentry_dsn
    },

    {
      name  = "SENTRY_DSN_ARN"
      value = aws_secretsmanager_secret.sentry_dsn.arn
    },

    {
      name  = "BATCH_EXECUTION_ROLE_ARN"
      value = aws_iam_role.data_refinery_batch_execution.arn
    },
    {
      name  = "S3_BUCKET_NAME"
      value = aws_s3_bucket.data_refinery_bucket.id
    },
    {
      name  = "S3_COMPENDIA_BUCKET_NAME"
      value = aws_s3_bucket.data_refinery_compendia_bucket.id
    },
    {
      name  = "S3_RESULTS_BUCKET_NAME"
      value = aws_s3_bucket.data_refinery_results_bucket.id
    },
    {
      name  = "S3_TRANSCRIPTOME_INDEX_BUCKET_NAME"
      value = aws_s3_bucket.data_refinery_transcriptome_index_bucket.id
    },
    {
      name  = "S3_QN_TARGET_BUCKET_NAME"
      value = aws_s3_bucket.data_refinery_qn_target_bucket.id
    },
    {
      name  = "LOCAL_ROOT_DIR"
      value = var.local_root_dir
    },
    {
      name  = "DOCKERHUB_REPO"
      value = var.dockerhub_repo
    },
    {
      name  = "DOWNLOADERS_DOCKER_IMAGE"
      value = var.downloaders_docker_image
    },
    {
      name  = "TRANSCRIPTOME_DOCKER_IMAGE"
      value = var.transcriptome_docker_image
    },
    {
      name  = "SALMON_DOCKER_IMAGE"
      value = var.salmon_docker_image
    },
    {
      name  = "AFFYMETRIX_DOCKER_IMAGE"
      value = var.affymetrix_docker_image
    },
    {
      name  = "ILLUMINA_DOCKER_IMAGE"
      value = var.illumina_docker_image
    },
    {
      name  = "NO_OP_DOCKER_IMAGE"
      value = var.no_op_docker_image
    },
    {
      name  = "FOREMAN_DOCKER_IMAGE"
      value = var.foreman_docker_image
    },
    {
      name  = "SMASHER_DOCKER_IMAGE"
      value = var.smasher_docker_image
    },
    {
      name  = "API_DOCKER_IMAGE"
      value = var.api_docker_image
    },
    {
      name  = "MAX_JOBS_PER_NODE"
      value = var.max_jobs_per_node
    },
    {
      name  = "MAX_DOWNLOADER_JOBS_PER_NODE"
      value = var.max_downloader_jobs_per_node
    },
    {
      name  = "SLACK_WEBHOOK_URL"
      value = var.slack_webhook_url
    },
    {
      name  = "REFINEBIO_JOB_QUEUE_WORKERS_NAMES"
      value = module.batch.data_refinery_workers_queue_names
    },
    {
      name  = "REFINEBIO_JOB_QUEUE_SMASHER_NAME"
      value = module.batch.data_refinery_smasher_queue_name
    },
    {
      name  = "REFINEBIO_JOB_QUEUE_COMPENDIA_NAME"
      value = module.batch.data_refinery_compendia_queue_name
    },
    {
      name  = "REFINEBIO_JOB_QUEUE_ALL_NAMES"
      value = module.batch.data_refinery_all_queue_names
    },
  ]
}
