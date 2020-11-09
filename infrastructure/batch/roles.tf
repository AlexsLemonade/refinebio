# This sets specific roles for access to AWS services
# used by compute environments and batch queues


### Batch Role
resource "aws_iam_role" "data_refinery_batch_role" {
  name = "data-refinery-batch-service-role"
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


# I can't tell if ECS is optional or not. I'm hoping so but I'm leaving this in case its not.
# ### ECS Role
# resource "aws_iam_role" "data_refinery_ecs_role" {
#   name = "data-refinery-ecs-instance-role"
#   assume_role_policy = <<EOF
# {
#     "Version": "2012-10-17",
#     "Statement": [
#     {
#         "Action": "sts:AssumeRole",
#         "Effect": "Allow",
#         "Principal": {
#         "Service": "ec2.amazonaws.com"
#         }
#     }
#     ]
# }
# EOF
#   tags = var.default_tags
# }

# resource "aws_iam_role_policy_attachment" "ecs_ec2_container" {
#   role = aws_iam_role.data_refinery_ecs_role.name
#   policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
# }

# resource "aws_iam_role_policy_attachment" "ecs_rw_s3" {
#   role = aws_iam_role.data_refinery_ecs_role.name
#   policy_arn = aws_iam_policy.data_refinery_readwrite_S3.arn
# }

# resource "aws_iam_role_policy_attachment" "ecs_read_s3" {
#   role = aws_iam_role.data_refinery_ecs_role.name
#   policy_arn = aws_iam_policy.data_refinery_read_S3.arn
# }


# ### Spotfleet Role
# resource "aws_iam_role" "data_refinery_spotfleet_role" {
#   name = "data-refinery-spotfleet-role"
#   assume_role_policy = <<EOF
# {
#     "Version": "2012-10-17",
#     "Statement": [
#     {
#         "Action": "sts:AssumeRole",
#         "Effect": "Allow",
#         "Principal": {
#         "Service": "spotfleet.amazonaws.com"
#         }
#     }
#     ]
# }
# EOF
#   tags = var.default_tags
# }

# resource "aws_iam_role_policy_attachment" "data_refinery_spotfleet_tagging" {
#   role = aws_iam_role.data_refinery_spotfleet_role.name
#   policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
# }

# resource "aws_iam_role_policy_attachment" "data_refinery_spotfleet_autoscale" {
#   role = aws_iam_role.data_refinery_spotfleet_role.name
#   policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetAutoscaleRole"
# }
