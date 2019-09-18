# Variables
# Only lowercase alphanumeric characters and hyphens allowed!
# In CI, this can be set with TF_VAR_user, TF_VAR_stage, etc.

# Annoyingly, TF can't have computed variables (${env.USER})
# as default values for variables. So `TF_VAR_user=rjones tf plan`, etc.
# Also, don't use any non-alphanumeric characters here or RDS will whinge.
variable "user" {
  default = "myusername"
}

variable "stage" {
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
variable "nomad_server_instance_type" {
  default = "m5.xlarge"
}

variable "client_instance_type" {
  default = "m5.4xlarge"
}

variable "smasher_instance_type" {
  # 256 GiB Memory, compendia jobs needs 220 and smasher jobs need 28.
  default = "m5.16xlarge"
}

variable "spot_price" {
  default = "4.10"
}

variable "spot_fleet_capacity" {
  default = "0"
}

variable "max_clients" {
  default = "12"
}

variable "scale_up_threshold" {
  default = "1"
}

variable "scale_down_threshold" {
  default = "0"
}

variable "raven_dsn" {
  default = ""
}

variable "raven_dsn_api" {
  default = ""
}

# API
variable "api_docker_image" {
  default = "dr_api:latest"
}

variable "api_instance_type" {
  default = "t2.large"
}

variable "foreman_instance_type" {
  default = "m5.large"
}

variable "volume_size_in_gb" {
  default = "9000"
}

variable "max_downloader_jobs_per_node" {
  default = 8
}

variable "elasticsearch_port" {
  default = "80"
}

variable "engagementbot_webhook" {
}

# Output our production environment variables.
output "environment_variables" {
  value = [
    {name = "REGION"
      value = "${var.region}"},
    {name = "USER"
      value = "${var.user}"},
    {name = "STAGE"
      value = "${var.stage}"},
    {name = "AWS_ACCESS_KEY_ID_CLIENT"
      value = "${aws_iam_access_key.data_refinery_user_client_key.id}"},
    {name = "AWS_SECRET_ACCESS_KEY_CLIENT"
      value = "${aws_iam_access_key.data_refinery_user_client_key.secret}"},
    {name = "DJANGO_DEBUG"
      value = "${var.django_debug}"},
    {name = "DJANGO_SECRET_KEY"
      value = "${var.django_secret_key}"},
    {name = "DATABASE_NAME"
      value = "${aws_db_instance.postgres_db.name}"},
    {name = "DATABASE_HOST"
      value = "${aws_instance.pg_bouncer.private_ip}"},
    {name = "RDS_HOST"
      value = "${aws_db_instance.postgres_db.address}"},
    {name = "DATABASE_USER"
      value = "${var.database_user}"},
    {name = "DATABASE_PASSWORD"
      value = "${var.database_password}"},
    {name = "DATABASE_PORT"
      value = "${var.database_port}"},
    {name = "DATABASE_HIDDEN_PORT"
      value = "${var.database_hidden_port}"},
    {name = "DATABASE_TIMEOUT"
      value = "${var.database_timeout}"},
    {name = "ELASTICSEARCH_HOST"
      value = "${aws_elasticsearch_domain.es.endpoint}"},
    {name = "ELASTICSEARCH_PORT"
      value = "${var.elasticsearch_port}"},
    {name = "RUNNING_IN_CLOUD"
      value = "${var.running_in_cloud}"},
    {name = "LOG_LEVEL"
      value = "${var.log_level}"},
    {name = "USE_S3"
      value = "${var.use_s3}"},
    {name = "RAVEN_DSN"
      value = "${var.raven_dsn}"},
    {name = "RAVEN_DSN_API"
      value = "${var.raven_dsn_api}"},
    {name = "S3_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_bucket.id}"},
    {name = "S3_RESULTS_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_results_bucket.id}"},
    {name = "S3_TRANSCRIPTOME_INDEX_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_transcriptome_index_bucket.id}"},
    {name = "S3_QN_TARGET_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_qn_target_bucket.id}"},
    {name = "S3_COMPENDIA_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_compendia_bucket.id}"},
    {name = "LOCAL_ROOT_DIR"
      value = "${var.local_root_dir}"},
    {name = "DOCKERHUB_REPO"
      value = "${var.dockerhub_repo}"},
    {name = "DOWNLOADERS_DOCKER_IMAGE"
      value = "${var.downloaders_docker_image}"},
    {name = "TRANSCRIPTOME_DOCKER_IMAGE"
      value = "${var.transcriptome_docker_image}"},
    {name = "SALMON_DOCKER_IMAGE"
      value = "${var.salmon_docker_image}"},
    {name = "AFFYMETRIX_DOCKER_IMAGE"
      value = "${var.affymetrix_docker_image}"},
    {name = "ILLUMINA_DOCKER_IMAGE"
      value = "${var.illumina_docker_image}"},
    {name = "NO_OP_DOCKER_IMAGE"
      value = "${var.no_op_docker_image}"},
    {name = "FOREMAN_DOCKER_IMAGE"
      value = "${var.foreman_docker_image}"},
    {name = "SMASHER_DOCKER_IMAGE"
      value = "${var.smasher_docker_image}"},
    {name = "API_DOCKER_IMAGE"
      value = "${var.api_docker_image}"},
    {name = "NOMAD_HOST"
      value = "${aws_instance.nomad_server_1.private_ip}"},
    {name = "NOMAD_PUBLIC_HOST"
      value = "${aws_instance.nomad_server_1.public_ip}"},
    {name = "NOMAD_PORT"
      value = "4646"},
    {name = "MAX_CLIENTS"
      value = "${var.max_clients}"},
    {name = "MAX_DOWNLOADER_JOBS_PER_NODE"
      value = "${var.max_downloader_jobs_per_node}"},
    {name = "ENGAGEMENTBOT_WEBHOOK"
      value = "${var.engagementbot_webhook}"}
  ]
}
