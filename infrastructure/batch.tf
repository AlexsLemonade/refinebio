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

  data_refinery_worker_user_data = data.template_file.nomad_client_script_smusher.rendered
  data_refinery_worker_ami = var.worker_ami

  user = var.user
  stage = var.stage
  region = var.region
}
