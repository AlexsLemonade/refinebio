# This file creates a nextflow user and adds it to the nextflow group.
# This user can be used for job submissions, but jobs can also be submitted as individual users.

resource "aws_iam_user" "data_refinery_user" {
  name = "data-refinery-batch-user"
  tags = var.default_tags
}

resource "aws_iam_access_key" "data_refinery_key" {
  user = aws_iam_user.data_refinery_user.name
}

resource "aws_iam_user_group_membership" "data_refinery_batch" {
  user = aws_iam_user.data_refinery_user.name
  groups = [
    aws_iam_group.data_refinery_group.name
  ]
}

resource "aws_iam_policy" "batch_policy" {
  name = "data-refinery-batch-policy-${var.user}-${var.stage}"
  description = "A policy that gives access to submit Batch jobs."
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
         "Effect":"Allow",
         "Action":[
            "batch:SubmitJob"
         ],
         "Resource":"arn:aws:batch:::*"
      }
    ]
}
EOF
}

resource "aws_iam_user_policy_attachment" "attach_batch_policy" {
  user = aws_iam_user.data_refinery_user.name
  policy_arn = aws_iam_policy.batch_policy.arn
}
