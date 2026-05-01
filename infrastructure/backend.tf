terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.25.0"
    }
    local = {
      source  = "hashicorp/local"
      version = "~> 2.4.0"
    }
  }
  backend "s3" {
    bucket = "refinebio-tfstate-deploy-${var.stage}"
    encrypt = true
    key = "terraform-${var.user}.tfstate"
    region = "us-east-1"
  }
}
