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
        "Service": ["ec2.amazonaws.com"]
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

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

}

resource "aws_iam_policy_attachment" "fleet_role" {
  name = "EC2SpotFleetRole"
  roles = [aws_iam_role.data_refinery_spot_fleet.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
}

resource "aws_iam_instance_profile" "data_refinery_instance_profile" {
  name = "data-refinery-instance-profile-${var.user}-${var.stage}"
  role = aws_iam_role.data_refinery_instance.name
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
          "Resource": [
            "arn:aws:s3:::${aws_s3_bucket.data_refinery_bucket.bucket}/*",
            "arn:aws:s3:::${aws_s3_bucket.data_refinery_results_bucket.bucket}/*",
            "arn:aws:s3:::${aws_s3_bucket.data_refinery_transcriptome_index_bucket.bucket}/*",
            "arn:aws:s3:::${aws_s3_bucket.data_refinery_qn_target_bucket.bucket}/*",
            "arn:aws:s3:::${aws_s3_bucket.data_refinery_compendia_bucket.bucket}/*"
          ]
      }
   ]
}
EOF

}

resource "aws_iam_role_policy_attachment" "s3" {
  role = aws_iam_role.data_refinery_instance.name
  policy_arn = aws_iam_policy.s3_access_policy.arn
}

resource "aws_iam_policy" "ec2_access_policy" {
  name = "data-refinery-ec2-access-policy-${var.user}-${var.stage}"
  description = "Allows EC2 Permissions."
  count = var.full_stack ? 1 : 0

  # We can't iterate instances from the fleet, so allow attaching to any instance,
  # but restrict which volumes can be attached.
  policy = <<EOF
{
   "Version":"2012-10-17",
   "Statement":[
      {
         "Effect":"Allow",
         "Action": [
            "ec2:DescribeVolumes"
          ],
          "Resource": [
            "*"
          ]
      },
      {
         "Effect":"Allow",
         "Action": [
            "ec2:AttachVolume"
          ],
          "Resource": [
            "${aws_ebs_volume.data_refinery_ebs_smasher[0].arn}",
            "${aws_instance.smasher_instance[0].arn}"
          ]
      },
      {
        "Effect": "Allow",
        "Action": [
          "sts:DecodeAuthorizationMessage"
        ],
        "Resource": [
          "*"
        ]
      }
   ]
}
EOF

}

resource "aws_iam_role_policy_attachment" "ec2" {
  role = aws_iam_role.data_refinery_instance.name
  policy_arn = aws_iam_policy.ec2_access_policy[0].arn
  count = var.full_stack ? 1 : 0
}

resource "aws_iam_policy" "cloudwatch_policy" {
  name = "data-refinery-cloudwatch-policy-${var.user}-${var.stage}"
  description = "Allows Cloudwatch Permissions."

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/iam-identity-based-access-control-cwl.html

  # Log streams are created dynamically by Nomad, so we give permission to the entire group
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
                "logs:DescribeLogStreams",
                "cloudwatch:GetMetricStatistics",
                "cloudwatch:ListMetrics",
                "cloudwatch:PutMetricAlarm",
                "cloudwatch:PutMetricData",
                "cloudwatch:SetAlarmState"
            ],
            "Resource": [
              "arn:aws:logs:${var.region}:${data.aws_caller_identity.current.account_id}:log-group:${aws_cloudwatch_log_group.data_refinery_log_group.name}:*"
            ]
        }
    ]
}
EOF

}

resource "aws_iam_role_policy_attachment" "cloudwatch" {
  role = aws_iam_role.data_refinery_instance.name
  policy_arn = aws_iam_policy.cloudwatch_policy.arn
}

resource "aws_s3_bucket_policy" "cloudtrail_access_policy" {
  bucket = aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.id
  policy = <<POLICY
{
  "Statement": [
    {
      "Action": "s3:GetBucketAcl",
      "Effect": "Allow",
      "Principal": {
        "Service": "cloudtrail.amazonaws.com"
      },
      "Resource": "${aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.arn}",
      "Sid": "AWSCloudTrailAclCheck20150319"
    },
    {
      "Action": "s3:PutObject",
      "Effect": "Allow",
      "Principal": {
        "Service": "cloudtrail.amazonaws.com"
      },
      "Resource": "${aws_s3_bucket.data_refinery_cloudtrail_logs_bucket.arn}/*",
      "Sid": "AWSCloudTrailWrite20150319"
    }
  ],
  "Version": "2012-10-17"
}
POLICY

}
