# This file creates the compute environments used by the queues
# The default environment is a 2 vCPU spot cluster

# Create an spot instance environment with up to 4 vcpus
# the AMI used is described in setup-log.md
resource "aws_batch_compute_environment" "data_refinery_spot" {
  compute_environment_name = "data-refinery-spot-compute-${var.user}-${var.stage}"
  compute_resources {
    instance_role = aws_iam_instance_profile.ecs_instance_profile.arn
    instance_type = [
      "m5.12xlarge", "r5.12xlarge"
    ]
    # allocation_strategy = "BEST_FIT_PROGRESSIVE"
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100
    max_vcpus = 48
    min_vcpus = 0
    # standard launch template
    launch_template {
      launch_template_id = aws_launch_template.data_refinery_launch_template.id
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

  service_role = aws_iam_role.batch_service_role.arn
  type = "MANAGED"
  depends_on = [aws_iam_role_policy_attachment.batch_service_role]
}
