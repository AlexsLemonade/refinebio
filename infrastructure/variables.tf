# Variables
# Only lowercase alphanumeric characters and hyphens allowed!
# In CI, this can be set with TF_VAR_user, TF_VAR_stage, etc.

# Annoyingly, TF can't have computed variables (${env.USER})
# as default values for variables. So `TF_VAR_user=rjones tf plan`, etc.
variable "user" {
  default = "my-user-name"
}

variable "stage" {
  default = "dev"
}

variable "region" {
  default = "us-east-1"
}
