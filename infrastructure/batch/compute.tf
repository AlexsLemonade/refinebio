# This file creates the compute environments used by the default and priority queues
# The default environment is a 100 vCPU spot cluster
# Priority environment is a 20 vCPU on demand cluster

# This may or may not prove necessary once more debugging can be done.
# resource "aws_iam_instance_profile" "data_refinery_ecs_instance_role" {
#   name = "data-refinery-ecs-instance-role"
#   role = aws_iam_role.data_refinery_ecs_role.name
# }

# Create an spot instance environment with up to 256 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "data_refinery_spot" {
  compute_environment_name = "data-refinery-spot-compute"
  compute_resources {
    instance_role = var.data_refinery_instance_profile.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100
    max_vcpus = 256
    min_vcpus = 0
    # standard launch template
    launch_template {
      launch_template_id = aws_launch_template.data_refinery_lt_standard.id
    }
    # ec2_key_pair = aws_key_pair.data_refinery_keypair.key_name
    security_group_ids = [
      aws_security_group.data_refinery_security.id,
    ]
    subnets = [
      aws_subnet.data_refinery_subnet.id,
    ]
    type = "SPOT"
    tags = merge(
      var.default_tags,
      {
        parent = "data-refinery-spot-compute"
      }
    )
  }

  service_role = aws_iam_role.data_refinery_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.data_refinery_batch_role]
}

# Create an spot instance0 environment with up to 32 vcpus with large disks

resource "aws_batch_compute_environment" "data_refinery_spot_bigdisk" {
  compute_environment_name = "data-refinery-spot-compute-bigdisk"
  compute_resources {
    instance_role = var.data_refinery_instance_profile.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100
    max_vcpus = 32
    min_vcpus = 0
    # large disk launch template
    launch_template {
      launch_template_id = aws_launch_template.data_refinery_lt_bigdisk.id
    }
    # ec2_key_pair = aws_key_pair.data_refinery_keypair.key_name
    security_group_ids = [
      aws_security_group.data_refinery_security.id,
    ]
    subnets = [
      aws_subnet.data_refinery_subnet.id,
    ]
    type = "SPOT"
    tags = merge(
      var.default_tags,
      {
        parent = "data-refinery-spot-compute"
      }
    )
  }

  service_role = aws_iam_role.data_refinery_batch_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.data_refinery_batch_role]
}

# # Create an ondemand environment with up to 32 vcpus
# # the AMI used is described in setup-log.md
# resource "aws_batch_compute_environment" "data_refinery_ondemand" {
#   compute_environment_name = "data-refinery-ondemand-compute"
#   compute_resources {
#     instance_role = var.data_refinery_instance_profile.arn
#     instance_type = [
#       "optimal",
#     ]
#     allocation_strategy = "BEST_FIT"
#     max_vcpus = 32
#     min_vcpus = 0
#     # standard launch template
#     launch_template {
#       launch_template_id = aws_launch_template.data_refinery_lt_standard.id
#     }
#     # ec2_key_pair = aws_key_pair.data_refinery_keypair.key_name
#     security_group_ids = [
#       aws_security_group.data_refinery_security.id,
#     ]
#     subnets = [
#       aws_subnet.data_refinery_subnet.id,
#     ]
#     type = "EC2"
#     tags = merge(
#       var.default_tags,
#       {
#         parent = "data-refinery-ondemand-compute"
#       }
#     )
#   }

#   service_role = aws_iam_role.data_refinery_batch_role.arn
#   type         = "MANAGED"
#   depends_on   = [aws_iam_role_policy_attachment.data_refinery_batch_role]

# }
