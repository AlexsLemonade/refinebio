resource "aws_iam_access_key" "data_refinery_user_worker_key" {
  user    = "${aws_iam_user.data_refinery_user_worker.name}"
}

resource "aws_iam_access_key" "data_refinery_user_foreman_key" {
  user    = "${aws_iam_user.data_refinery_user_foreman.name}"
}

resource "aws_iam_access_key" "data-refinery-deployer-access-key" {
  user = "${aws_iam_user.data-refinery-deployer.name}"
}

resource "aws_iam_user" "data_refinery_user_worker" {
  name = "data-refinery-user-worker-${var.user}-${var.stage}"
}

resource "aws_iam_user" "data_refinery_user_foreman" {
  name = "data-refinery-user-foreman-${var.user}-${var.stage}"
}

resource "aws_iam_user" "data-refinery-deployer" {
  name = "data-refinery-deployer-${var.user}-${var.stage}"
}


# XXX: TODO: Lock these down!!!!
resource "aws_iam_user_policy" "data_refinery_user_worker_policy" {
  name = "data-refinery-user-worker-key-${var.user}-${var.stage}"
  user = "${aws_iam_user.data_refinery_user_worker.name}"

  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
          		  "logs:PutLogEvents",
                "logs:DescribeLogStreams"
            ],
            "Resource": "arn:aws:logs:*:*:*"
        },
        {
            "Effect": "Allow",
            "Action": [
                "s3:*"
            ],
            "Resource": "arn:aws:s3:::*"
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
        },
	{
	    "Effect":"Allow",
	    "Action":[
		"ses:SendEmail",
		"ses:SendRawEmail"
	    ],
	    "Resource":"*"
	}
    ]
}
EOF
}

resource "aws_iam_user_policy" "data_refinery_user_foreman_policy" {
  name = "data-refinery-user-foreman-key-${var.user}-${var.stage}"
  user = "${aws_iam_user.data_refinery_user_foreman.name}"

  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
		            "logs:PutLogEvents",
		            "logs:DescribeLogStreams"
            ],
            "Resource":
                "arn:aws:logs:*:*:*"
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

data "aws_iam_policy_document" "data-refinery-deployment" {
  statement {
    actions = [
      "s3:ListObjects",
      "s3:GetObject",
      "s3:PutObject",
      "s3:DeleteObject",
      "s3:PutObjectAcl",
    ]

    resources = [
      "arn:aws:s3:::${aws_s3_bucket.data-refinery-static.id}/*",
    ]
  }

  statement {
    actions = [
      "s3:ListBucket"
    ]

    resources = [
      "arn:aws:s3:::${aws_s3_bucket.data-refinery-static.id}",
    ]
  }
}

resource "aws_iam_user_policy" "data-refinery-deployer" {
  name = "data-refinery-deployer-${var.user}-${var.stage}"
  user = "${aws_iam_user.data-refinery-deployer.name}"
  policy = "${data.aws_iam_policy_document.data-refinery-deployment.json}"
}
