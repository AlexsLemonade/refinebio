resource "aws_iam_access_key" "data_refinery_user_client_key" {
  user    = "${aws_iam_user.data_refinery_user_client.name}"
}

resource "aws_iam_access_key" "data-refinery-deployer-access-key" {
  user = "${aws_iam_user.data-refinery-deployer.name}"
}

resource "aws_iam_user" "data_refinery_user_client" {
  name = "data-refinery-user-client-${var.user}-${var.stage}"
}

resource "aws_iam_user" "data-refinery-deployer" {
  name = "data-refinery-deployer-${var.user}-${var.stage}"
}


# XXX: TODO: Lock these down!!!!
resource "aws_iam_user_policy" "data_refinery_user_client_policy" {
  name = "data-refinery-user-client-key-${var.user}-${var.stage}"
  user = "${aws_iam_user.data_refinery_user_client.name}"

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
            "Resource": [
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_surveyor.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_processor.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_downloader.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_foreman.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_api.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_api_nginx_access.name}",
              "arn:aws:logs:${var.region}:${aws_cloudwatch_log_group.data_refinery_log_group.name}:${aws_cloudwatch_log_stream.log_stream_api_nginx_error.name}"
            ]
        },
        {
            "Effect": "Allow",
            "Action": [
              "s3:AbortMultipartUpload",
              "s3:DeleteObject",
              "s3:DeleteObjectTagging",
              "s3:DeleteObjectVersion",
              "s3:GetObject",
              "s3:GetObjectAcl",
              "s3:GetObjectTagging",
              "s3:GetObjectVersion",
              "s3:GetObjectVersionAcl",
              "s3:ListMultipartUploadParts",
              "s3:PutObject",
              "s3:PutObjectAcl"
            ],
            "Resource": [
              "arn:aws:s3:::${aws_s3_bucket.data_refinery_bucket.bucket}/*",
              "arn:aws:s3:::${aws_s3_bucket.data_refinery_results_bucket.bucket}/*",
              "arn:aws:s3:::${aws_s3_bucket.data_refinery_transcriptome_index_bucket.bucket}/*",
              "arn:aws:s3:::${aws_s3_bucket.data_refinery_compendia_bucket.bucket}/*"
            ]
        },
        {
            "Effect": "Allow",
            "Action":[
              "SES:SendEmail",
              "SES:SendRawEmail"
            ],
            "Resource": "arn:aws:ses:${var.region}:${data.aws_caller_identity.current.account_id}:identity/refine.bio"
        },
        {
            "Effect": "Allow",
            "Action":[
              "ec2:DescribeVolumes",
              "ec2:AttachVolume"
            ],
            "Resource": [
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 0)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 1)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 2)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 3)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 4)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 5)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 6)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 7)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 8)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 9)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 10)}",
              "arn:aws:ec2:${var.region}:${data.aws_caller_identity.current.account_id}:volume/${element(aws_ebs_volume.data_refinery_ebs.*.id, 11)}",
            ]
        },
        {
            "Effect": "Allow",
            "Action":[
              "es:DescribeElasticsearchDomain",
              "es:DescribeElasticsearchDomainConfig",
              "es:DescribeElasticsearchDomains",
              "es:DescribeElasticsearchInstanceTypeLimits",
              "es:ListDomainNames",
              "es:ListElasticsearchInstanceTypeDetails",
              "es:ListElasticsearchInstanceTypes",
              "es:ListElasticsearchVersions",
              "es:ListTags",
              "es:ESHttpDelete",
              "es:ESHttpGet",
              "es:ESHttpHead",
              "es:ESHttpPost",
              "es:ESHttpPut"
            ],
            "Resource": "arn:aws:es:${var.region}:${data.aws_caller_identity.current.account_id}:domain/${aws_elasticsearch_domain.es.domain_name}"
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
