# Variables
# Only lowercase alphanumeric characters and hyphens allowed!
# In CI, this can be set with TF_VAR_user, TF_VAR_stage, etc.

# Annoyingly, TF can't have computed variables (${env.USER})
# as default values for variables. So `TF_VAR_user=rjones tf plan`, etc.
# Also, don't use any non-alphanumeric characters here or RDS will whinge.
variable "user" {
  default = "rjones"
}

variable "stage" {
  default = "dev"
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
  # Not sure what's going on with this here password, but it doesn't seem to work.
  default = "drpostgrespassword"
}

variable "django_secret_key" {
  default = "NtG1bxZU115GThwrLuAJe0PhTVN9hJ4P"
}
variable "django_debug" {
  default = "False"
}
variable "database_port" {
  default = "5432"
}

variable "database_timeout" {
  default = "30"
}

variable "running_in_cloud" {
  default = "True"
}

# This is a placeholder until there is a production image ready.
variable "workers_docker_image" {
  default = "wkurt/dr_worker:3"
}
variable "foreman_docker_image" {
  default = "wkurt/dr_foreman:4"
}
variable "use_s3" {
  default = "True"
}
variable "s3_bucket_name" {
  default = "data-refinery"
}
variable "local_root_dir" {
  default = "/home/user/data_store"
}
variable "raw_prefix" {
  default = "raw"
}
variable "temp_prefix" {
  default = "temp"
}
variable "processed_prefix" {
  default = "processed"
}

# Output our production environment variables. These should be in
# parity with the env files such as workers/environments/dev.
output "environment_variables" {
  value = [
    {name = "REGION"
      value = "${var.region}"},
    {name = "USER"
      value = "${var.user}"},
    {name = "STAGE"
      value = "${var.stage}"},
    {name = "AWS_ACCESS_KEY_ID"
      value = "${aws_iam_access_key.data_refinery_user_worker_key.id}"},
    {name = "AWS_SECRET_ACCESS_KEY"
      value = "${aws_iam_access_key.data_refinery_user_worker_key.secret}"},
    {name = "DJANGO_DEBUG"
      value = "${var.django_debug}"},
    {name = "DJANGO_SECRET_KEY"
      value = "${var.django_secret_key}"},
    {name = "DATABASE_NAME"
      value = "${aws_db_instance.postgres_db.name}"},
    {name = "DATABASE_HOST"
      value = "${aws_db_instance.postgres_db.address}"},
    {name = "DATABASE_USER"
      value = "${var.database_user}"},
    {name = "DATABASE_PASSWORD"
      value = "${var.database_password}"},
    {name = "DATABASE_PORT"
      value = "${var.database_port}"},
    {name = "DATABASE_TIMEOUT"
      value = "${var.database_timeout}"},
    {name = "RUNNING_IN_CLOUD"
      value = "${var.running_in_cloud}"},
    {name = "USE_S3"
      value = "${var.use_s3}"},
    {name = "S3_BUCKET_NAME"
      value = "${aws_s3_bucket.data_refinery_bucket.id}"},
    {name = "LOCAL_ROOT_DIR"
      value = "${var.local_root_dir}"},
    {name = "RAW_PREFIX"
      value = "${var.raw_prefix}"},
    {name = "TEMP_PREFIX"
      value = "${var.temp_prefix}"},
    {name = "PROCESSED_PREFIX"
      value = "${var.processed_prefix}"},
    {name = "WORKERS_DOCKER_IMAGE"
      value = "${var.workers_docker_image}"},
    {name = "FOREMAN_DOCKER_IMAGE"
      value = "${var.foreman_docker_image}"},
    {name = "NOMAD_HOST"
      value = "${aws_instance.nomad_server_1.private_ip}"},
    {name = "NOMAD_PUBLIC_HOST"
      value = "${aws_instance.nomad_server_1.public_ip}"},
    {name = "NOMAD_PORT"
      value = "4646"}
  ]
}
