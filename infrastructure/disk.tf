# Elastic File System and EFS Mount Points, and S3

##
# EFS
##

resource "aws_efs_file_system" "data_refinery_efs" {
  creation_token = "data-refinery-efs-${var.user}-${var.stage}"

  tags {
    Name = "data-refinery-efs-${var.user}-${var.stage}"
  }
}

resource "aws_efs_mount_target" "data_refinery_efs_1a" {
  file_system_id = "${aws_efs_file_system.data_refinery_efs.id}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  security_groups = ["${aws_security_group.data_refinery_worker.id}"]
}

resource "aws_efs_mount_target" "data_refinery_efs_1b" {
  file_system_id = "${aws_efs_file_system.data_refinery_efs.id}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  security_groups = ["${aws_security_group.data_refinery_worker.id}"]
}

##
# EBS
## 

resource "aws_ebs_volume" "data_refinery_ebs" {
  count = "${var.max_clients}"
  availability_zone = "${var.region}a"
  size = 1600 # 16TB
  type = "st1" # Throughput Optimized HDD
  tags {
    Name        = "data-refinery-ebs-${count.index}-${var.user}-${var.stage}"
    Environment = "${var.stage}"
    Index       = "${count.index"}"
  }
}

##
# S3
##

resource "aws_s3_bucket" "data_refinery_bucket" {
  bucket = "data-refinery-s3-${var.user}-${var.stage}"
  acl    = "private"

  tags {
    Name        = "data-refinery-s3-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }
}

resource "aws_s3_bucket" "data_refinery_results_bucket" {
  bucket = "data-refinery-s3-results-${var.user}-${var.stage}"
  acl    = "public-read"

  tags {
    Name        = "data-refinery-s3-results-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }

  lifecycle_rule {
    id = "auto-delete-after-1-days-${var.user}-${var.stage}"
    prefix = ""
    enabled = true
    abort_incomplete_multipart_upload_days = 1

    expiration {
      days = 1
      expired_object_delete_marker = true
    }

    noncurrent_version_expiration {
      days = 1
    }
  }
}

resource "aws_s3_bucket" "data-refinery-static" {
  bucket = "${var.static_bucket_prefix == "dev" ? var.user : var.static_bucket_prefix}${var.static_bucket_root}"

  cors_rule {
    allowed_origins = ["*"]
    allowed_methods = ["GET"]
    max_age_seconds = 3000
    allowed_headers = ["Authorization"]
  }

  tags {
    Name = "data-refinery-static-site-${var.user}-${var.stage}"
  }

  website {
    index_document = "index.html"
  }
}

resource "aws_s3_bucket" "data_refinery_transcriptome_index_bucket" {
  bucket = "data-refinery-s3-transcriptome-index-${var.user}-${var.stage}"
  acl    = "private"

  tags {
    Name        = "data-refinery-s3-transcriptome-index-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }
}

resource "aws_s3_bucket" "data-refinery-static-access-logs" {
  bucket = "data-refinery-static-access-logs-${var.user}-${var.stage}"

  tags {
    Name = "data-refinery-static-access-logs-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }

  lifecycle_rule {
    id = "auto-delete-after-1-days-${var.user}-${var.stage}"
    prefix = ""
    enabled = true
    abort_incomplete_multipart_upload_days = 7

    expiration {
      days = 7
      expired_object_delete_marker = true
    }

    noncurrent_version_expiration {
      days = 7
    }
  }
}
