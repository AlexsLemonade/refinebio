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

variable "region" {
  default = "us-east-1"
}

variable "database_user" {
  default = "drpostgresuser"
}

variable "database_password" {
  default = "drpostgrespassword"
}