resource "aws_security_group" "ami_template_instance" {
  name = "data-refinery-ami-template-instance"
  description = "data-refinery-ami-template-instance"
  vpc_id = data.aws_vpc.ccdl_dev_vpc.id

  tags {
    Name = "data-refinery-ami-template-instance"
  }
}

# Allow ssh connections
resource "aws_security_group_rule" "ssh" {
  type = "ingress"
  from_port = 22
  to_port = 22
  protocol = "tcp"
  cidr_blocks = ["0.0.0.0/0"]
  security_group_id = aws_security_group.ami_template_instance.id
}

# Allow all outbound traffic
resource "aws_security_group_rule" "outbound" {
  type = "egress"
  from_port = 0
  to_port = 65535
  protocol = "all"
  cidr_blocks = ["0.0.0.0/0"]
  ipv6_cidr_blocks = ["::/0"]
  security_group_id = aws_security_group.ami_template_instance.id
}
