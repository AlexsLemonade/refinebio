# EBS Volumes and S3

##
# EBS
##

resource "aws_ebs_volume" "data_refinery_ebs" {
  count = "${var.max_clients}"
  availability_zone = "${var.region}a"
  size = "${var.volume_size_in_gb}" # 1600 = 1.6TB
  type = "st1" # Throughput Optimized HDD
  tags {
    Name        = "data-refinery-ebs-${count.index}-${var.user}-${var.stage}"
    Environment = "${var.stage}"
    Index       = "${count.index}"
    User        = "${var.user}"
    Stage       = "${var.stage}"
    IsBig       = "True"
  }
}

##
# S3
##

resource "aws_s3_bucket" "data_refinery_bucket" {
  bucket = "data-refinery-s3-${var.user}-${var.stage}"
  acl    = "private"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

  tags {
    Name        = "data-refinery-s3-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }
}

resource "aws_s3_bucket" "data_refinery_results_bucket" {
  bucket = "data-refinery-s3-results-${var.user}-${var.stage}"
  acl    = "private"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

  tags {
    Name        = "data-refinery-s3-results-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }

  lifecycle_rule {
    id = "auto-delete-after-7-days-${var.user}-${var.stage}"
    prefix = ""
    enabled = true
    abort_incomplete_multipart_upload_days = 1

    expiration {
      days = 7
      expired_object_delete_marker = true
    }

    noncurrent_version_expiration {
      days = 1
    }
  }
}

resource "aws_s3_bucket" "data-refinery-static" {
  bucket = "${var.static_bucket_prefix == "dev" ? var.user : var.static_bucket_prefix}${var.static_bucket_root}"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

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
  acl    = "public-read"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

  tags {
    Name        = "data-refinery-s3-transcriptome-index-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }
}

resource "aws_s3_bucket" "data_refinery_qn_target_bucket" {
  bucket = "data-refinery-s3-qn-target-${var.user}-${var.stage}"
  acl    = "public-read"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

  tags {
    Name        = "data-refinery-s3-qn-target-${var.user}-${var.stage}"
    Environment = "${var.stage}"
  }
}

resource "aws_s3_bucket" "data_refinery_compendia_bucket" {
  bucket = "data-refinery-s3-compendia-${var.user}-${var.stage}"
  acl    = "private"
  force_destroy = "${var.static_bucket_prefix == "dev" ? true : false}"

  tags {
    Name        = "data-refinery-s3-compendia-${var.user}-${var.stage}"
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
    id = "auto-delete-after-7-days-${var.user}-${var.stage}"
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
