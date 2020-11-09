resource "aws_launch_template" "data_refinery_lt_standard" {
  name = "data-refinery-launchtemplate-standard"
  tags = var.default_tags
  # ccdl-data-refinery-base-v2.0 image
  image_id = "ami-01c08d4de548477df"
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 64
      encrypted = true
      delete_on_termination = true
    }
  }
}

resource "aws_launch_template" "data_refinery_lt_bigdisk" {
  name = "data-refinery-launchtemplate-bigdisk"
  tags = var.default_tags
  # ccdl-data-refinery-base-v2.0 image
  image_id = "ami-01c08d4de548477df"
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 500
      encrypted = true
      delete_on_termination = true
    }
  }
}
