terraform {
  backend "s3" {
    bucket = "refinebio-tfstate-deploy-${var.stage}"
    encrypt = true
    key = "terraform-${var.user}.tfstate"
    region = "us-east-1"
  }
}
