resource "aws_ecs_cluster" "data_refinery" {
  name = "data-refinery"
}

resource "aws_iam_role" "ecs_service_role" {
  name = "data-refinery-ecs-service"

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
  deployment_minimum_healthy_percent = 50
  deployment_maximum_percent = 200
  iam_role = "${aws_iam_role.ecs_service_role.name}"
  depends_on = ["aws_iam_role_policy.ecs_service"]

  load_balancer {
    elb_name = "${aws_elb.data_refinery_worker.name}"
    container_name = "data-refinery-worker"
    container_port = 8000
  }
}
