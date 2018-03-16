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
resource "aws_cloudwatch_log_stream" "log_stream_nomad_server_1" {
  name           = "log-stream-nomad-server-1-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_nomad_server_2" {
  name           = "log-stream-nomad-server-2-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}
resource "aws_cloudwatch_log_stream" "log_stream_nomad_server_3" {
  name           = "log-stream-nomad-server-3-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# XXX: Once we're using spot instnaces we should have them create the
# streams themselves so that we'll have streams for the dynamically
# created instances. However, in the meantime explicitly creating the
# stream via Terraform is better because it will be able to clean them
# up too.
resource "aws_cloudwatch_log_stream" "log_stream_nomad_client" {
  name           = "log-stream-nomad-client-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

# Nomad / Docker
resource "aws_cloudwatch_log_stream" "log_stream_nomad_docker_surveyor" {
  name           = "log-stream-nomad-docker-surveyor-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_nomad_docker_processor" {
  name           = "log-stream-nomad-docker-processor-${var.user}-${var.stage}"
  log_group_name = "${aws_cloudwatch_log_group.data_refinery_log_group.name}"
}

resource "aws_cloudwatch_log_stream" "log_stream_nomad_docker_downloader" {
  name           = "log-stream-nomad-docker-downloader-${var.user}-${var.stage}"
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
