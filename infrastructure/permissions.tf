# The configuration contained in this file specifies AWS IAM roles and
# permissions.

resource "aws_iam_role" "data_refinery_instance" {
  name = "data-refinery-instance-${var.user}-${var.stage}"

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
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF
}

resource "aws_iam_instance_profile" "data_refinery_instance_profile" {
  name  = "data-refinery-instance-profile-${var.user}-${var.stage}"
  role = "${aws_iam_role.data_refinery_instance.name}"
}

resource "aws_iam_policy" "s3_access_policy" {
  name = "data-refinery-s3-access-policy-${var.user}-${var.stage}"
  description = "Allows S3 Permissions."

  # Policy text based off of:
  # http://docs.aws.amazon.com/AmazonS3/latest/dev/example-bucket-policies.html
  policy = <<EOF
{
   "Version":"2012-10-17",
   "Statement":[
      {
         "Effect":"Allow",
         "Action":[
            "s3:ListAllMyBuckets"
         ],
         "Resource":"arn:aws:s3:::*"
      },
      {
         "Effect":"Allow",
         "Action":[
            "s3:ListBucket",
            "s3:GetBucketLocation"
         ],
         "Resource":"arn:aws:s3:::data-refinery"
      },
      {
         "Effect":"Allow",
         "Action":[
            "s3:PutObject",
            "s3:GetObject",
            "s3:DeleteObject"
         ],
         "Resource":"arn:aws:s3:::data-refinery/*"
      }
   ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "s3" {
  role = "${aws_iam_role.data_refinery_instance.name}"
  policy_arn = "${aws_iam_policy.s3_access_policy.arn}"
}

resource "aws_iam_policy" "ec2_access_policy" {
  name = "data-refinery-ec2-access-policy-${var.user}-${var.stage}"
  description = "Allows EC2 Permissions."

  # Policy text based off of:
  # http://docs.aws.amazon.com/AmazonS3/latest/dev/example-bucket-policies.html
  policy = <<EOF
{
   "Version":"2012-10-17",
   "Statement":[
      {
         "Effect":"Allow",
         "Action":[
            "ec2:*"
         ],
         "Resource": "*"
      }
   ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "ec2" {
  role = "${aws_iam_role.data_refinery_instance.name}"
  policy_arn = "${aws_iam_policy.ec2_access_policy.arn}"
}

resource "aws_iam_policy" "cloudwatch_policy" {
  name = "data-refinery-cloudwatch-policy-${var.user}-${var.stage}"
  description = "Allows Cloudwatch Permissions."


  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/iam-identity-based-access-control-cwl.html
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "logs:CreateLogGroup",
                "logs:CreateLogStream",
                "logs:PutLogEvents",
                "logs:DescribeLogStreams"
            ],
            "Resource": [
                "arn:aws:logs:*:*:*"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
                "cloudwatch:GetMetricStatistics",
                "cloudwatch:ListMetrics",
                "cloudwatch:PutMetricAlarm",
                "cloudwatch:PutMetricData",
                "cloudwatch:SetAlarmState"
            ],
            "Resource":
                "*"
        }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "cloudwatch" {
  role = "${aws_iam_role.data_refinery_instance.name}"
  policy_arn = "${aws_iam_policy.cloudwatch_policy.arn}"
}
