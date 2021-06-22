# AWS Batch setup
provider "aws" {
  profile = "default"
  region = "us-east-1"
}

resource "aws_batch_job_queue" "data_refinery_workers_queues" {
  count = var.num_workers

  name = "data-refinery-batch-workers-queue-${var.user}-${var.stage}-${count.index}"
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_workers[count.index].arn,
  ]

  tags = var.default_tags
}

output "data_refinery_workers_queue_names" {
  value = join(",", local.worker_queue_names)
}

locals {
  worker_queue_names = [for q in aws_batch_job_queue.data_refinery_workers_queues : q.name]
}

resource "aws_batch_job_queue" "data_refinery_smasher_queue" {
  name = "data-refinery-batch-smasher-queue-${var.user}-${var.stage}"
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_smasher.arn,
  ]

  tags = var.default_tags
}

output "data_refinery_smasher_queue_name" {
  value = aws_batch_job_queue.data_refinery_smasher_queue.name
}

resource "aws_batch_job_queue" "data_refinery_compendia_queue" {
  name = "data-refinery-batch-compendia-queue-${var.user}-${var.stage}"
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.data_refinery_compendia.arn,
  ]

  tags = var.default_tags
}

output "data_refinery_compendia_queue_name" {
  value = aws_batch_job_queue.data_refinery_compendia_queue.name
}

locals {
  all_queue_names = concat(
    local.worker_queue_names,
    [aws_batch_job_queue.data_refinery_smasher_queue.name,
     aws_batch_job_queue.data_refinery_compendia_queue.name]
  )
}

output "data_refinery_all_queue_names" {
  value = join(",", local.all_queue_names)
}
