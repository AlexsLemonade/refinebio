provider "aws" {
  region = "us-east-1"
}

resource "aws_s3_bucket" "refinebio-tfstate" {
  bucket = "refinebio-tfstate-${var.user}-${var.stage}"
  acl    = "private"

  versioning {
    enabled = true
  }
}

output "terraform_state_s3_bucket" {
  value = "${aws_s3_bucket.refinebio-tfstate.id}"
}
