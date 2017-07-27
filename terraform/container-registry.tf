# The configuration contained in this file creates an ECS service,
# cluster, and task definitions for The Data Refinery. It also creates
# the IAM permissions needed for ECS to work.

resource "aws_ecs_cluster" "data_refinery" {
  name = "data-refinery"
}

# The following IAM role/policies allow ECS to register/deregister EC2
# instances with the ELB. More information on why these are needed can
# be found at:
# http://docs.aws.amazon.com/AmazonECS/latest/developerguide/service_IAM_role.html
resource "aws_iam_role" "ecs_service_role" {
  name = "data-refinery-ecs-service"

  # Similar policy text can be found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/service_IAM_role.html
  # I am unsure of where this exact text came from, but it's nearly identical.
  assume_role_policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": "sts:AssumeRole",
      "Principal": {
        "Service": "ecs.amazonaws.com"
      },
      "Effect": "Allow",
      "Sid": ""
    }
  ]
}
EOF
}

resource "aws_iam_role_policy" "ecs_service" {
  name = "data-refinery-ecs-service-policy"
  role = "${aws_iam_role.ecs_service_role.name}"

  # Policy text found at:
  # http://docs.aws.amazon.com/AmazonECS/latest/developerguide/service_IAM_role.html
  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "ec2:AuthorizeSecurityGroupIngress",
                "ec2:Describe*",
                "elasticloadbalancing:DeregisterInstancesFromLoadBalancer",
                "elasticloadbalancing:DeregisterTargets",
                "elasticloadbalancing:Describe*",
                "elasticloadbalancing:RegisterInstancesWithLoadBalancer",
                "elasticloadbalancing:RegisterTargets"
            ],
            "Resource": "*"
        }
    ]
}
EOF
}

resource "aws_ecs_task_definition" "data_refinery_worker" {
  family = "data-refinery-worker"
  container_definitions = "${file("task-definitions/worker.json.secret")}"
}

resource "aws_ecs_service" "data_refinery_worker" {
  name = "data-refinery-worker"
  cluster = "${aws_ecs_cluster.data_refinery.id}"
  task_definition = "${aws_ecs_task_definition.data_refinery_worker.arn}"
  desired_count  = 2
  deployment_minimum_healthy_percent = 0
  deployment_maximum_percent = 200
  depends_on = ["aws_iam_role_policy.ecs_service"]

  lifecycle {
    ignore_changes = ["task_definition"]
  }
}
