resource "aws_iam_access_key" "data_refinery_user_worker_key" {
  user    = "${aws_iam_user.data_refinery_user_worker.name}"
}

resource "aws_iam_access_key" "data_refinery_user_foreman_key" {
  user    = "${aws_iam_user.data_refinery_user_foreman.name}"
}

resource "aws_iam_user" "data_refinery_user_worker" {
  name = "data-refinery-user-worker-${var.user}-${var.stage}"
}

resource "aws_iam_user" "data_refinery_user_foreman" {
  name = "data-refinery-user-foreman-${var.user}-${var.stage}"
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
        }
    ]
}
EOF
}
