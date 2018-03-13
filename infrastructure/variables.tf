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
}

variable "database_user" {
  default = "drpostgresuser"
}

variable "database_password" {
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

# Can be hardcoded into application
variable "nomad_host" {
  default = "nomad"
}

# This is a placeholder until there is a production image ready.
variable "workers_docker_image" {
  default = "miserlou/dr_worker:2"
}
variable "foreman_docker_image" {
  default = "miserlou/dr_foreman:3"
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
