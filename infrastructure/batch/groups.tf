# This file sets the access settings for the data-refinery-batch group
# Currently includes read acess to all s3, and write access to select buckets

resource "aws_iam_group" "data_refinery_group" {
  name = "data-refinery-batch-${var.user}-${var.stage}"
}

# Batch access
resource "aws_iam_group_policy_attachment" "batch_access" {
  group = aws_iam_group.data_refinery_group.name
  policy_arn = "arn:aws:iam::aws:policy/AWSBatchFullAccess"
}
