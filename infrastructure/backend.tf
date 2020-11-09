terraform {
  backend "s3" {
    # Terraform will prompt the user for the other keys.
    region = "us-east-1"
  }
}

data "terraform_remote_state" "network" {
  backend = "s3"
  config = {
    bucket = "refinebio-tfstate-deploy-${var.stage}"
    key = "terraform-${var.user}.tfstate"
    region = "us-east-1"
    encrypt = true
  }
}
