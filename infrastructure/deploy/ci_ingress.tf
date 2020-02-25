resource "aws_security_group_rule" "data_refinery_ci_nomad" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  cidr_blocks = ["${var.host_ip}/32"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_security_group_rule" "data_refinery_ci_postgres" {
  type = "ingress"
  from_port = "${var.database_port}"
  to_port = "${var.database_port}"
  protocol = "tcp"
  cidr_blocks = ["${var.host_ip}/32"]
  security_group_id = "${aws_security_group.data_refinery_pg.id}"
}
