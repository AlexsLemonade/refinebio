# VPC and network security settings

resource "aws_security_group" "data_refinery_security" {
  name = "data-refinery-security-group"
  vpc_id = aws_vpc.data_refinery_vpc.id
  tags = var.default_tags

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  # ingress {
  #   description = "SSH from anywhere."
  #   from_port   = 22
  #   to_port     = 22
  #   protocol    = "tcp"
  #   cidr_blocks = ["0.0.0.0/0"]
  # }
}

# resource "aws_key_pair" "data_refinery_keypair" {
#   key_name = "data-refinery-key"
#   public_key = "PUT_YOUR_PUBLIC_KEY_HERE"
# }

resource "aws_vpc" "data_refinery_vpc" {
  cidr_block = "10.1.0.0/16"
  tags = merge(
    {
      Name = "data-refinery-vpc"
    },
    var.default_tags
  )
  enable_dns_hostnames = true
}

resource "aws_internet_gateway" "data_refinery_gateway" {
  vpc_id = aws_vpc.data_refinery_vpc.id

  tags = merge(
    {
      Name = "data-refinery-gateway"
    },
    var.default_tags
  )

}

resource "aws_subnet" "data_refinery_subnet" {
  vpc_id = aws_vpc.data_refinery_vpc.id
  cidr_block = "10.1.1.0/24"
  tags = merge(
    {
      Name = "data-refinery-subnet"
    },
    var.default_tags
  )
  map_public_ip_on_launch = true
}

resource "aws_route_table" "data_refinery_route_table" {
  vpc_id = aws_vpc.data_refinery_vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.data_refinery_gateway.id
  }

  tags = merge(
    {
      Name = "data-refinery-route-table"
    },
    var.default_tags
  )
}

resource "aws_route_table_association" "data_refinery_route_table_association" {
  subnet_id = aws_subnet.data_refinery_subnet.id
  route_table_id = aws_route_table.data_refinery_route_table.id
}
