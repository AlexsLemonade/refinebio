# The configuration contained in this file specifies AWS resources
# related to networking.

provider "aws" {
  region = var.region
}

resource "aws_vpc" "data_refinery_vpc" {
  cidr_block = "10.0.0.0/16"
  enable_dns_support = true
  enable_dns_hostnames = true

  tags = merge(
    var.default_tags,
    {
      Name = "data-refinery-${var.user}-${var.stage}"
    }
  )
}

resource "aws_subnet" "data_refinery_1a" {
  availability_zone = "${var.region}a"
  cidr_block = "10.0.0.0/17"
  vpc_id = aws_vpc.data_refinery_vpc.id
  map_public_ip_on_launch = true

  tags = merge(
    var.default_tags,
    {
      Name = "data-refinery-1a-${var.user}-${var.stage}"
    }
  )
}

resource "aws_subnet" "data_refinery_1b" {
  availability_zone = "${var.region}b"
  cidr_block = "10.0.128.0/17"
  vpc_id = aws_vpc.data_refinery_vpc.id

  # Unsure if this should be set to true
  map_public_ip_on_launch = true

  tags = merge(
    var.default_tags,
    {
      Name = "data-refinery-1b-${var.user}-${var.stage}"
    }
  )
}

resource "aws_internet_gateway" "data_refinery" {
  vpc_id = aws_vpc.data_refinery_vpc.id

  tags = merge(
    var.default_tags,
    {
      Name = "data-refinery-${var.user}-${var.stage}"
    }
  )
}

# Note: this is a insecure practice long term, however it's
# necessary to access it from lab machines.
resource "aws_route_table" "data_refinery" {
  vpc_id = aws_vpc.data_refinery_vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.data_refinery.id
  }

  tags = merge(
    var.default_tags,
    {
      Name = "data-refinery-${var.user}-${var.stage}"
    }
  )
}

resource "aws_route_table_association" "data_refinery_1a" {
  subnet_id = aws_subnet.data_refinery_1a.id
  route_table_id = aws_route_table.data_refinery.id
}

resource "aws_route_table_association" "data_refinery_1b" {
  subnet_id = aws_subnet.data_refinery_1b.id
  route_table_id = aws_route_table.data_refinery.id
}

resource "aws_db_subnet_group" "data_refinery" {
  name = "data-refinery-${var.user}-${var.stage}"
  subnet_ids = [aws_subnet.data_refinery_1a.id, aws_subnet.data_refinery_1b.id]

  tags = merge(
    var.default_tags,
    {
      Name = "Data Refinery DB Subnet ${var.user}-${var.stage}"
    }
  )
}

# Get the API a static IP address.
resource "aws_eip" "data_refinery_api_ip" {
  vpc = true

  tags = merge(
    var.default_tags,
    {
      Name = "Data Refinery API Elastic IP ${var.user}-${var.stage}"
    }
  )
}

# As per https://aws.amazon.com/elasticloadbalancing/details/:
#
# You can select the appropriate load balancer based on your
# application needs. If you need flexible application management, we
# recommend that you use an Application Load Balancer. If extreme
# performance and static IP is needed for your application, we
# recommend that you use a Network Load Balancer. If you have an
# existing application that was built within the EC2-Classic network,
# then you should use a Classic Load Balancer.
#
# it appears an Network Load Balancer would be best for us because we
# need a static IP address to point DNS to.
resource "aws_lb" "data_refinery_api_load_balancer" {
  # Extra short because there is a 32 char limit on this name
  name = "DR-api-${var.user}-${var.stage}"
  internal = false
  load_balancer_type = "network"

  # Only one subnet is allowed and the API lives in 1a.
  subnet_mapping {
    subnet_id = aws_subnet.data_refinery_1a.id
    allocation_id = aws_eip.data_refinery_api_ip.id
  }

  tags = var.default_tags
}

resource "aws_lb_target_group" "api-http" {
  name = "dr-api-${var.user}-${var.stage}-http"
  port = 80
  protocol = "TCP"
  vpc_id = aws_vpc.data_refinery_vpc.id

  tags = var.default_tags
}

resource "aws_lb_listener" "api-http" {
  load_balancer_arn = aws_lb.data_refinery_api_load_balancer.arn
  protocol = "TCP"
  port = 80

  default_action {
    target_group_arn = aws_lb_target_group.api-http.arn
    type = "forward"
  }
}

resource "aws_lb_target_group_attachment" "api-http" {
  target_group_arn = aws_lb_target_group.api-http.arn
  target_id = aws_instance.api_server_1.id
  port = 80
}

resource "aws_lb_target_group" "api-https" {
  name = "dr-api-${var.user}-${var.stage}-https"
  port = 443
  protocol = "TCP"
  vpc_id = aws_vpc.data_refinery_vpc.id

  tags = var.default_tags
}

resource "aws_lb_listener" "api-https" {
  load_balancer_arn = aws_lb.data_refinery_api_load_balancer.arn
  protocol = "TCP"
  port = 443

  default_action {
    target_group_arn = aws_lb_target_group.api-https.arn
    type = "forward"
  }
}

resource "aws_lb_target_group_attachment" "api-https" {
  target_group_arn = aws_lb_target_group.api-https.arn
  target_id = aws_instance.api_server_1.id
  port = 443
}
