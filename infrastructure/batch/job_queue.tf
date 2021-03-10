# AWS Batch setup
provider "aws" {
  profile = "default"
  region = "us-east-1"
}

variable "default_tags" {
  description = "Default resource tags"
  type = map(string)
  default = {
    purpose = "data-refinery-batch"
  }

}

resource "aws_batch_job_queue" "data_refinery_default_queue" {
  name = "data-refinery-batch-default-queue-${var.user}-${var.stage}"
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_spot.arn,
  ]
}

output "data_refinery_default_queue_name" {
  value = aws_batch_job_queue.data_refinery_default_queue.name
}
