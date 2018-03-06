# The configuration contained in this file specifies AWS resources
# related to networking.

provider "aws" {
  region = "${var.region}"
}

resource "aws_vpc" "data_refinery_vpc" {
  cidr_block = "10.0.0.0/16"
  enable_dns_support = true
  enable_dns_hostnames = true

  tags {
    Name = "data-refinery-${var.user}-${var.stage}"
  }
}

resource "aws_subnet" "data_refinery_1a" {
  availability_zone = "${var.region}a"
  cidr_block = "10.0.0.0/17"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"
  map_public_ip_on_launch = true

  tags {
    Name = "data-refinery-1a-${var.user}-${var.stage}"
  }
}

resource "aws_subnet" "data_refinery_1b" {
  availability_zone = "${var.region}b"
  cidr_block = "10.0.128.0/17"
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"
  # Unsure if this should be set to true
  map_public_ip_on_launch = true

  tags {
    Name = "data-refinery-1b-${var.user}-${var.stage}"
  }
}

resource "aws_internet_gateway" "data_refinery" {
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  tags = {
    Name = "data-refinery-${var.user}-${var.stage}"
  }
}

# Note: this is a insecure practice long term, however it's
# necessary to access it from lab machines.
resource "aws_route_table" "data_refinery" {
  vpc_id = "${aws_vpc.data_refinery_vpc.id}"

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = "${aws_internet_gateway.data_refinery.id}"
  }

  tags {
    Name = "data-refinery-${var.user}-${var.stage}"
  }
}


resource "aws_route_table_association" "data_refinery_1a" {
  subnet_id      = "${aws_subnet.data_refinery_1a.id}"
  route_table_id = "${aws_route_table.data_refinery.id}"
}

resource "aws_route_table_association" "data_refinery_1b" {
  subnet_id      = "${aws_subnet.data_refinery_1b.id}"
  route_table_id = "${aws_route_table.data_refinery.id}"
}

resource "aws_db_subnet_group" "data_refinery" {
  name = "data-refinery-${var.user}-${var.stage}"
  subnet_ids = ["${aws_subnet.data_refinery_1a.id}", "${aws_subnet.data_refinery_1b.id}"]

  tags {
    Name = "Data Refinery DB Subnet ${var.user}-${var.stage}"
  }
}
