module "batch" {
  source = "./batch"

  # VPC
  data_refinery_vpc = aws_vpc.data_refinery_vpc
  data_refinery_subnet = aws_subnet.data_refinery_1a

  # S3 Buckets
  data_refinery_bucket = aws_s3_bucket.data_refinery_bucket
  data_refinery_results_bucket = aws_s3_bucket.data_refinery_results_bucket
  data_refinery_transcriptome_index_bucket = aws_s3_bucket.data_refinery_transcriptome_index_bucket
  data_refinery_qn_target_bucket = aws_s3_bucket.data_refinery_qn_target_bucket
  data_refinery_compendia_bucket = aws_s3_bucket.data_refinery_compendia_bucket

  data_refinery_spot_fleet_role = aws_iam_role.data_refinery_spot_fleet
  data_refinery_keypair = aws_key_pair.data_refinery
  data_refinery_worker_security_group = aws_security_group.data_refinery_worker

  data_refinery_worker_user_data = templatefile(
    "workers-configuration/workers-instance-user-data.tpl.sh",
    {
      database_host     = aws_instance.pg_bouncer.private_ip
      database_name     = aws_db_instance.postgres_db.name
      database_password = var.database_password
      database_port     = var.database_port
      database_user     = var.database_user
      region            = var.region
      stage             = var.stage
      user              = var.user
    }
  )
  data_refinery_worker_ami = var.worker_ami

  user = var.user
  stage = var.stage
  region = var.region
  default_tags = var.default_tags
  num_workers = var.num_workers
  use_on_demand_instances = var.batch_use_on_demand_instances
}
