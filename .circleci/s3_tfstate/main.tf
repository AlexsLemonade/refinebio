provider "aws" {
  region = "us-east-1"
}


resource "aws_s3_bucket" "refinebio-tfstate" {
  bucket = "refinebio-tfstate"
  acl    = "private"

  versioning {
    enabled = true
  }
}
