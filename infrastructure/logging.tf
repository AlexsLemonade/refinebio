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

##
# Metrics and Alarms
##

# Turns out we don't need to actually _make_ the metric - we can just push to it
# without creating it via TF and it JustWorks^tm - we just have to make sure
# that we match to this value in our script.

# We need one metric for scaling up and another for scaling down.
# If the queue length is larger than [threshold], on average,
# for more than [evaluation_periods] of [period] seconds, then
# fire the alarm action to the scale up command, which adds a server.

# For the down alarm, do the opposite.

# resource "aws_cloudwatch_metric_alarm" "nomad_queue_length_alarm_up" {
#     alarm_name = "nomad-queue-length-alarm-up-${var.user}-${var.stage}"
#     comparison_operator = "GreaterThanOrEqualToThreshold"
#     evaluation_periods = "2"
#     metric_name = "NomadQueueLength"
#     namespace = "${var.user}-${var.stage}"
#     period = "120"
#     statistic = "Average"
#     threshold = "${var.scale_up_threshold}"
#     alarm_description = "The queue is too long - we need more workers!"
#     alarm_actions = [
#         "${aws_autoscaling_policy.clients_scale_up.arn}"
#     ]

# }

# resource "aws_cloudwatch_metric_alarm" "nomad_queue_length_alarm_down" {
#     alarm_name = "nomad-queue-length-alarm-down-${var.user}-${var.stage}"
#     comparison_operator = "LessThanOrEqualToThreshold"
#     evaluation_periods = "2"
#     metric_name = "NomadQueueLength"
#     namespace = "${var.user}-${var.stage}"
#     period = "120"
#     statistic = "Average"
#     threshold = "${var.scale_down_threshold}"
#     alarm_description = "The queue is too short - we need less workers!"
#     alarm_actions = [
#         "${aws_autoscaling_policy.clients_scale_down.arn}"
#     ]

# }
