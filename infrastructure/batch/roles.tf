# This sets specific roles for access to AWS services
# used by compute environments and batch queues


### Batch Role
resource "aws_iam_role" "data_refinery_batch_role" {
  name = "data-refinery-batch-service-role-${var.user}-${var.stage}"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "batch.amazonaws.com"
        }
    }
    ]
}
EOF
  tags = var.default_tags
}

resource "aws_iam_role_policy_attachment" "data_refinery_batch_role" {
  role = aws_iam_role.data_refinery_batch_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}
