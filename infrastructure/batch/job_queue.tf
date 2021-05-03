# AWS Batch setup
provider "aws" {
  profile = "default"
  region = "us-east-1"
}

resource "aws_batch_job_queue" "data_refinery_default_queue" {
  name = "data-refinery-batch-default-queue-${var.user}-${var.stage}"
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_spot.arn,
  ]

  tags = var.default_tags
}

output "data_refinery_default_queue_name" {
  value = aws_batch_job_queue.data_refinery_default_queue.name
}
