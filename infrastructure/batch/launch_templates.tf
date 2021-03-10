resource "aws_launch_template" "data_refinery_lt_standard" {
  name = "data-refinery-launchtemplate-standard-${var.user}-${var.stage}"
  tags = var.default_tags
  # ccdl-data-refinery-base-v2.0 image
  image_id = var.data_refinery_worker_ami
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 64
      encrypted = true
      delete_on_termination = true
    }
  }

  user_data = base64encode(var.data_refinery_worker_user_data)
}

resource "aws_launch_template" "data_refinery_lt_bigdisk" {
  name = "data-refinery-launchtemplate-bigdisk-${var.user}-${var.stage}"
  tags = var.default_tags
  # ccdl-data-refinery-base-v2.0 image
  image_id = var.data_refinery_worker_ami
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size = 500
      encrypted = true
      delete_on_termination = true
    }
  }

  user_data = base64encode(var.data_refinery_worker_user_data)
}
