# CloudTrail

resource "aws_cloudtrail" "data_refinery_s3_cloudtrail" {
  name = "data-refinery-s3-cloudtrail-${var.user}-${var.stage}"
  depends_on = ["aws_s3_bucket.data_refinery_cloudtrail_logs_bucket"]
  s3_bucket_name = "${aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id}"
  include_global_service_events = false
  event_selector {
    read_write_type = "ReadOnly"
    include_management_events = false

    data_resource {
      type = "AWS::S3::Object"

      # Make sure to append a trailing '/' to your ARN if you want
      # to monitor all objects in a bucket.
      # ref https://www.terraform.io/docs/providers/aws/r/cloudtrail.html#logging-individual-s3-bucket-events
      values = [
        "${aws_s3_bucket.data_refinery_compendia_bucket.arn}/",
        "${aws_s3_bucket.data_refinery_transcriptome_index_bucket.arn}/",
        "${aws_s3_bucket.data_refinery_qn_target_bucket.arn}/",
      ]
    }
  }
}

# CloudWatch Log Groups and Streams

##
# Groups
##

# This is the group. All of the streams belong to this.
resource "aws_cloudwatch_log_group" "data_refinery_log_group" {
  name = "data-refinery-log-group-${var.user}-${var.stage}"

  tags {
    Name = "data-refinery-log-group-${var.user}-${var.stage}"
  }
}

##
# Streams
##

# Nomad / Docker
resource "aws_cloudwatch_log_stream" "log_stream_surveyor" {
  name           = "log-stream-surveyor-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_processor" {
  name           = "log-stream-processor-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_downloader" {
  name           = "log-stream-downloader-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# Foreman
resource "aws_cloudwatch_log_stream" "log_stream_foreman" {
  name           = "log-stream-foreman-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# API
resource "aws_cloudwatch_log_stream" "log_stream_api" {
  name           = "log-stream-api-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_api_nginx_access" {
  name           = "log-stream-api-nginx-access-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_api_nginx_error" {
  name           = "log-stream-api-nginx-error-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# Cloudtrail Logs

# Must start with `/aws/events` in order to connect to a cloudwatch_event_target.
resource "aws_cloudwatch_log_group" "compendia_object_metrics_log_group" {
  name = "/aws/events/data-refinery-compendia-log-group-${var.user}-${var.stage}"
}
