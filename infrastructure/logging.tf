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

# Nomad
resource "aws_cloudwatch_log_stream" "log_stream_nomad_server" {
  name           = "log-stream-nomad-server-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_nomad_client" {
  name           = "log-stream-nomad-client-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# Foreman
resource "aws_cloudwatch_log_stream" "log_stream_foreman_django" {
  name           = "log-stream-foreman-django-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_foreman_docker" {
  name           = "log-stream-foreman-docker-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# Worker
resource "aws_cloudwatch_log_stream" "log_stream_worker_django" {
  name           = "log-stream-worker-django-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_worker_docker" {
  name           = "log-stream-worker-docker-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}
