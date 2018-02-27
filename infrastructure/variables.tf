# Variables
# Only lowercase alphanumeric characters and hyphens allowed!
# In CI, this can be set with TF_VAR_user, TF_VAR_stage, etc.

variable "user" {
  default = "my-user-name"
}

variable "stage" {
  default = "dev"
}