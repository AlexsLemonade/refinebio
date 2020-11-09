# AWS Batch setup
provider "aws" {
  profile = "default"
  region  = "us-east-1"
}

variable "default_tags" {
  description = "Default resource tags"
  type        = map(string)
  default     = {
    purpose = "data-refinery-batch"
    config = "https://github.com/AlexsLemonade/alsf-scpca/tree/jashapiro/terraform-batch/aws"
  }

}

resource "aws_batch_job_queue" "data_refinery_default_queue" {
  name     = "data-refinery-batch-default-queue"
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_spot.arn,
  ]
}

resource "aws_batch_job_queue" "data_refinery_bigdisk_queue" {
  name     = "data-refinery-batch-bigdisk-queue"
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_spot_bigdisk.arn,
  ]
}


# resource "aws_batch_job_queue" "data_refinery_priority_queue" {
#   name     = "data-refinery-batch-priority-queue"
#   state    = "ENABLED"
#   priority = 1
#   compute_environments = [
#     aws_batch_compute_environment.data_refinery_ondemand.arn,
#   ]
# }
