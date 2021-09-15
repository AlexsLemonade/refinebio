# This file creates the compute environments used by the queues
# The default environment is a 2 vCPU spot cluster

locals {
  on_demand_environment_type = "EC2"
  on_demand_allocation_strategy = "BEST_FIT"
  compute_environment_type = var.use_on_demand_instances == true ? local.on_demand_environment_type : "SPOT"
  allocation_strategy = var.use_on_demand_instances == true ? local.on_demand_allocation_strategy : "SPOT_CAPACITY_OPTIMIZED"
}

resource "aws_batch_compute_environment" "data_refinery_workers" {
  count = var.num_workers

  compute_environment_name = "data-refinery-workers-compute-${var.user}-${var.stage}-${count.index}"
  compute_resources {
    instance_role = aws_iam_instance_profile.ecs_instance_profile.arn
    instance_type = [
      # Large enough to be worth our time, small enough to be affordable.
      "m5.12xlarge", "m5.16xlarge", "r5.12xlarge", "r5.8xlarge"
    ]

    type = local.compute_environment_type
    allocation_strategy = local.allocation_strategy
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100

    max_vcpus = 32
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
    tags = merge(
      var.default_tags,
      {
        parent = "data-refinery-workers-compute"
      }
    )
  }

  service_role = aws_iam_role.batch_service_role.arn
  type = "MANAGED"
  depends_on = [aws_iam_role_policy_attachment.batch_service_role]

  tags = var.default_tags
}

# A separate queue for smasher jobs so they don't get stuck behind a
# large backlog.
resource "aws_batch_compute_environment" "data_refinery_smasher" {
  compute_environment_name = "data-refinery-smasher-compute-${var.user}-${var.stage}"
  compute_resources {
    instance_role = aws_iam_instance_profile.ecs_instance_profile.arn
    instance_type = [
      # Just large enough to run one job at a time. We haven't had a
      # backlog yet.
      "m5.2xlarge", "r5.xlarge"
    ]

    type = local.compute_environment_type
    allocation_strategy = local.allocation_strategy
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100

    max_vcpus = 32
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
    tags = merge(
      var.default_tags,
      {
        parent = "data-refinery-smasher-compute"
      }
    )
  }

  service_role = aws_iam_role.batch_service_role.arn
  type = "MANAGED"
  depends_on = [aws_iam_role_policy_attachment.batch_service_role]

  tags = var.default_tags
}

# Large instances for crunching hefty compendia, namely human and mouse.
resource "aws_batch_compute_environment" "data_refinery_compendia" {
  compute_environment_name = "data-refinery-compendia-compute-${var.user}-${var.stage}"
  compute_resources {
    instance_role = aws_iam_instance_profile.ecs_instance_profile.arn
    instance_type = [
      # 976 GB RAM
      "x1.16xlarge"
    ]

    # Compendia are long running, so we don't want to use spot instances.
    type = local.on_demand_environment_type
    allocation_strategy = local.on_demand_allocation_strategy
    spot_iam_fleet_role = var.data_refinery_spot_fleet_role.arn
    bid_percentage = 100

    max_vcpus = 32
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
    tags = merge(
      var.default_tags,
      {
        parent = "data-refinery-compendia-compute"
      }
    )
  }

  service_role = aws_iam_role.batch_service_role.arn
  type = "MANAGED"
  depends_on = [aws_iam_role_policy_attachment.batch_service_role]

  tags = var.default_tags
}
