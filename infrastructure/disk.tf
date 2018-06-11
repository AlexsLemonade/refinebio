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
  acl    = "public"

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
