resource "aws_s3_bucket" "data_refinery_bucket" {
  bucket        = "data-refinery-s3-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_bucket" {
  bucket = aws_s3_bucket.data_refinery_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_bucket" {
  bucket = aws_s3_bucket.data_refinery_bucket.id

  block_public_acls   = true
  block_public_policy = true
  ignore_public_acls  = true
}

resource "aws_s3_bucket_acl" "data_refinery_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_bucket,
    aws_s3_bucket_public_access_block.data_refinery_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_bucket.id
  acl    = "private"
}

resource "aws_s3_bucket" "data_refinery_results_bucket" {
  bucket        = "data-refinery-s3-results-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-results-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )

}

resource "aws_s3_bucket_ownership_controls" "data_refinery_results_bucket" {
  bucket = aws_s3_bucket.data_refinery_results_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_results_bucket" {
  bucket = aws_s3_bucket.data_refinery_results_bucket.id

  block_public_acls   = true
  block_public_policy = true
  ignore_public_acls  = true
}

resource "aws_s3_bucket_acl" "data_refinery_results_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_results_bucket,
    aws_s3_bucket_public_access_block.data_refinery_results_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_results_bucket.id
  acl    = "private"
}

resource "aws_s3_bucket_lifecycle_configuration" "data_refinery_results_bucket" {
  bucket = aws_s3_bucket.data_refinery_results_bucket.id

  rule {
    id     = "auto-delete-after-7-days-${var.user}-${var.stage}"
    status = "Enabled"

    abort_incomplete_multipart_upload {
      days_after_initiation = 1
    }

    expiration {
      days                         = 7
      expired_object_delete_marker = true
    }

    noncurrent_version_expiration {
      noncurrent_days = 1
    }
  }
}


resource "aws_s3_bucket" "data_refinery_transcriptome_index_bucket" {
  bucket        = "data-refinery-s3-transcriptome-index-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-transcriptome-index-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_transcriptome_index_bucket" {
  bucket = aws_s3_bucket.data_refinery_transcriptome_index_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_transcriptome_index_bucket" {
  bucket = aws_s3_bucket.data_refinery_transcriptome_index_bucket.id

  block_public_acls       = false
  block_public_policy     = false
  ignore_public_acls      = false
  restrict_public_buckets = false
}

resource "aws_s3_bucket_acl" "data_refinery_transcriptome_index_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_transcriptome_index_bucket,
    aws_s3_bucket_public_access_block.data_refinery_transcriptome_index_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_transcriptome_index_bucket.id
  acl    = "public-read"
}

resource "aws_s3_bucket" "data_refinery_qn_target_bucket" {
  bucket        = "data-refinery-s3-qn-target-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-qn-target-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_qn_target_bucket" {
  bucket = aws_s3_bucket.data_refinery_qn_target_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_qn_target_bucket" {
  bucket = aws_s3_bucket.data_refinery_qn_target_bucket.id

  block_public_acls       = false
  block_public_policy     = false
  ignore_public_acls      = false
  restrict_public_buckets = false
}

resource "aws_s3_bucket_acl" "data_refinery_qn_target_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_qn_target_bucket,
    aws_s3_bucket_public_access_block.data_refinery_qn_target_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_qn_target_bucket.id
  acl    = "public-read"
}

resource "aws_s3_bucket" "data_refinery_compendia_bucket" {
  bucket        = "data-refinery-s3-compendia-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-compendia-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_compendia_bucket" {
  bucket = aws_s3_bucket.data_refinery_compendia_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_compendia_bucket" {
  bucket = aws_s3_bucket.data_refinery_compendia_bucket.id

  block_public_acls   = true
  block_public_policy = true
  ignore_public_acls  = true
}

resource "aws_s3_bucket_acl" "data_refinery_compendia_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_compendia_bucket,
    aws_s3_bucket_public_access_block.data_refinery_compendia_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_compendia_bucket.id
  acl    = "private"
}


resource "aws_s3_bucket" "data_refinery_cloudtrail_logs_bucket" {
  bucket        = "data-refinery-s3-cloudtrail-logs-${var.user}-${var.stage}"
  force_destroy = var.static_bucket_prefix == "dev" ? true : false

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-s3-cloudtrail-logs-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_cloudtrail_logs_bucket" {
  bucket = aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_cloudtrail_logs_bucket" {
  bucket = aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id

  block_public_acls   = true
  block_public_policy = true
  ignore_public_acls  = true
}

resource "aws_s3_bucket_acl" "data_refinery_cloudtrail_logs_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_cloudtrail_logs_bucket,
    aws_s3_bucket_public_access_block.data_refinery_cloudtrail_logs_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id
  acl    = "private"
}


# Passing the name attribute as `EntireBucket` enables request metrics for the bucket.
# ref: https://www.terraform.io/docs/providers/aws/r/s3_bucket_metric.html
resource "aws_s3_bucket_metric" "compendia_bucket_metrics" {
  bucket = aws_s3_bucket.data_refinery_compendia_bucket.bucket
  name   = "EntireBucket"
}

resource "aws_cloudwatch_event_rule" "compendia_object_metrics" {
  name        = "data-refinery-compendia-object-metric-${var.user}-${var.stage}"
  description = "Download Compendia Events"

  event_pattern = <<PATTERN
{
  "source": [
    "aws.s3"
  ],
  "detail-type": [
    "AWS API Call via CloudTrail"
  ],
  "detail": {
    "eventSource": [
      "s3.amazonaws.com"
    ],
    "eventName": [
      "GetObject"
    ],
    "requestParameters": {
      "bucketName": [
        "${aws_s3_bucket.data_refinery_compendia_bucket.bucket}"
      ]
    }
  }
}
PATTERN

}

# arn needs to be stripped of trailing `:*`
# aws appends this when creating the event_target
resource "aws_cloudwatch_event_target" "compendia_object_metrics_target" {
  rule      = aws_cloudwatch_event_rule.compendia_object_metrics.name
  target_id = "compendia-object-logs-target-${var.user}-${var.stage}"
  arn       = aws_cloudwatch_log_group.compendia_object_metrics_log_group.arn
}

resource "aws_s3_bucket" "data_refinery_cert_bucket" {
  bucket        = "data-refinery-cert-${var.user}-${var.stage}"
  force_destroy = var.stage == "prod" ? false : true

  tags = merge(
    var.default_tags,
    {
      Name        = "data-refinery-cert-${var.user}-${var.stage}"
      Environment = var.stage
    }
  )
}

resource "aws_s3_bucket_ownership_controls" "data_refinery_cert_bucket" {
  bucket = aws_s3_bucket.data_refinery_cert_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_public_access_block" "data_refinery_cert_bucket" {
  bucket = aws_s3_bucket.data_refinery_cert_bucket.id

  block_public_acls   = true
  block_public_policy = true
  ignore_public_acls  = true
}

resource "aws_s3_bucket_acl" "data_refinery_cert_bucket" {
  depends_on = [
    aws_s3_bucket_ownership_controls.data_refinery_cert_bucket,
    aws_s3_bucket_public_access_block.data_refinery_cert_bucket,
  ]

  bucket = aws_s3_bucket.data_refinery_cert_bucket.id
  acl    = "private"
}

resource "aws_s3_bucket_lifecycle_configuration" "data_refinery_cert_bucket" {
  bucket = aws_s3_bucket.data_refinery_cert_bucket.id

  rule {
    id     = "auto-delete-after-30-days-${var.user}-${var.stage}"
    status = "Enabled"

    abort_incomplete_multipart_upload {
      days_after_initiation = 1
    }

    expiration {
      days                         = 30
      expired_object_delete_marker = true
    }

    noncurrent_version_expiration {
      noncurrent_days = 1
    }
  }
}
