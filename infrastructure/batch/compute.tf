# This file creates the compute environments used by the default and priority queues
# The default environment is a 4 vCPU spot cluster

# Create an spot instance environment with up to 4 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "data_refinery_spot" {
  compute_environment_name = "data-refinery-spot-compute-${var.user}-${var.stage}"
  compute_resources {
    instance_role = var.data_refinery_instance_profile.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100
    max_vcpus = 4
    min_vcpus = 0
    # standard launch template
    launch_template {
      launch_template_id = aws_launch_template.data_refinery_lt_standard.id
    }
    ec2_key_pair = var.data_refinery_keypair.key_name
    security_group_ids = [
      aws_security_group.data_refinery_security.id,
      var.data_refinery_worker_security_group.id,
    ]
    subnets = [
      var.data_refinery_subnet.id,
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
  compute_environment_name = "data-refinery-spot-compute-bigdisk-${var.user}-${var.stage}"
  compute_resources {
    instance_role = var.data_refinery_instance_profile.arn
    instance_type = [
      "optimal",
    ]
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100
    max_vcpus = 4
    min_vcpus = 0
    # large disk launch template
    launch_template {
      launch_template_id = aws_launch_template.data_refinery_lt_bigdisk.id
    }
    ec2_key_pair = var.data_refinery_keypair.key_name
    security_group_ids = [
      aws_security_group.data_refinery_security.id,
    ]
    subnets = [
      var.data_refinery_subnet.id,
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
