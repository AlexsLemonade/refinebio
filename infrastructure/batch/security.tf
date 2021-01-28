# VPC and network security settings

resource "aws_security_group" "data_refinery_security" {
  name = "data-refinery-security-group-${var.user}-${var.stage}"
  vpc_id = var.data_refinery_vpc.id
  tags = var.default_tags

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  ingress {
    description = "SSH from anywhere."
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }
}
