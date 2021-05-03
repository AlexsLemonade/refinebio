# The configuration contained in this file specifies AWS IAM roles and
# permissions.

resource "aws_iam_role" "data_refinery_worker" {
  name = "data-refinery-worker-${var.user}-${var.stage}"

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html
  # "Service": ["ec2.amazonaws.com", "ecs-tasks.amazonaws.com", "batch.amazonaws.com"]
  assume_role_policy = <<EOF
{
  "Version": "2008-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": ["ecs-tasks.amazonaws.com"]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

  tags = var.default_tags
}

resource "aws_iam_role" "data_refinery_foreman" {
  name = "data-refinery-foreman-${var.user}-${var.stage}"

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html
  assume_role_policy = <<EOF
{
  "Version": "2008-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": ["ec2.amazonaws.com"]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

  tags = var.default_tags
}

resource "aws_iam_role" "data_refinery_api" {
  name = "data-refinery-api-${var.user}-${var.stage}"

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/instance_IAM_role.html
  assume_role_policy = <<EOF
{
  "Version": "2008-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": ["ec2.amazonaws.com"]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

  tags = var.default_tags
}

resource "aws_iam_role" "data_refinery_spot_fleet" {
  name = "data-refinery-spot-fleet-${var.user}-${var.stage}"

  assume_role_policy = <<EOF
{
  "Version": "2008-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": ["spotfleet.amazonaws.com"]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

  tags = var.default_tags
}

resource "aws_iam_policy_attachment" "fleet_role" {
  name = "EC2SpotFleetRole"
  roles = [aws_iam_role.data_refinery_spot_fleet.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
}

resource "aws_iam_instance_profile" "data_refinery_worker" {
  name = "data-refinery-worker-profile-${var.user}-${var.stage}"
  role = aws_iam_role.data_refinery_worker.name
}

resource "aws_iam_instance_profile" "data_refinery_foreman" {
  name = "data-refinery-foreman-profile-${var.user}-${var.stage}"
  role = aws_iam_role.data_refinery_foreman.name
}

resource "aws_iam_instance_profile" "data_refinery_api" {
  name = "data-refinery-api-profile-${var.user}-${var.stage}"
  role = aws_iam_role.data_refinery_api.name
}

resource "aws_iam_policy" "s3_access_policy" {
  name = "data-refinery-s3-access-policy-${var.user}-${var.stage}"
  description = "Allows S3 Permissions."

  # Policy text based off of:
  # http://docs.aws.amazon.com/AmazonS3/latest/dev/example-bucket-policies.html
  policy = local.s3_access_policy

  tags = var.default_tags
}

resource "aws_iam_role_policy_attachment" "worker_s3" {
  role = aws_iam_role.data_refinery_worker.name
  policy_arn = aws_iam_policy.s3_access_policy.arn
}

resource "aws_iam_role_policy_attachment" "foreman_s3" {
  role = aws_iam_role.data_refinery_foreman.name
  policy_arn = aws_iam_policy.s3_access_policy.arn
}

resource "aws_iam_policy" "cloudwatch_policy" {
  name = "data-refinery-cloudwatch-policy-${var.user}-${var.stage}"
  description = "Allows Cloudwatch Permissions."

  policy = local.cloudwatch_logs_policy

  tags = var.default_tags
}

resource "aws_iam_role_policy_attachment" "foreman_cloudwatch" {
  role = aws_iam_role.data_refinery_foreman.name
  policy_arn = aws_iam_policy.cloudwatch_policy.arn
}

resource "aws_iam_role_policy_attachment" "api_cloudwatch" {
  role = aws_iam_role.data_refinery_api.name
  policy_arn = aws_iam_policy.cloudwatch_policy.arn
}

resource "aws_s3_bucket_policy" "cloudtrail_access_policy" {
  bucket = aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id
  policy = local.cloudtrail_access_policy
}

# Needed by Batch. Don't love that This has to be defined here instead of in Batch's module.
resource "aws_iam_role_policy_attachment" "ecs_ec2_container" {
  role = aws_iam_role.data_refinery_worker.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_policy" "batch_read_policy" {
  name = "data-refinery-batch-read-policy-${var.user}-${var.stage}"
  description = "Allows Batch Read Permissions."

  policy = local.batch_read_policy

  tags = var.default_tags
}

resource "aws_iam_policy" "batch_write_policy" {
  name = "data-refinery-batch-write-policy-${var.user}-${var.stage}"
  description = "Allows Batch Write Permissions."

  policy = local.batch_write_policy

  tags = var.default_tags
}

resource "aws_iam_role_policy_attachment" "worker_batch_read" {
  role = aws_iam_role.data_refinery_worker.name
  policy_arn = aws_iam_policy.batch_read_policy.arn
}

resource "aws_iam_role_policy_attachment" "worker_batch_write" {
  role = aws_iam_role.data_refinery_worker.name
  policy_arn = aws_iam_policy.batch_write_policy.arn
}

resource "aws_iam_role_policy_attachment" "foreman_batch_read" {
  role = aws_iam_role.data_refinery_foreman.name
  policy_arn = aws_iam_policy.batch_read_policy.arn
}

resource "aws_iam_role_policy_attachment" "foreman_batch_write" {
  role = aws_iam_role.data_refinery_foreman.name
  policy_arn = aws_iam_policy.batch_write_policy.arn
}

resource "aws_iam_role_policy_attachment" "api_batch_read" {
  role = aws_iam_role.data_refinery_api.name
  policy_arn = aws_iam_policy.batch_read_policy.arn
}

# The API needs to write to be able to submit smasher jobs
resource "aws_iam_role_policy_attachment" "api_batch_write" {
  role = aws_iam_role.data_refinery_api.name
  policy_arn = aws_iam_policy.batch_write_policy.arn
}
