resource "aws_launch_template" "data_refinery_worker" {
  name = "data-refinery-worker-${var.user}-${var.stage}"
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

  block_device_mappings {
    device_name = "/dev/sdf"
    ebs {
      volume_type = "st1"
      volume_size = 300
      encrypted = true
      delete_on_termination = true
    }
  }

  user_data = base64encode(var.data_refinery_worker_user_data)

  tag_specifications {
      resource_type = "instance"
      tags = var.default_tags
  }

  tag_specifications {
    resource_type = "volume"
    tags = var.default_tags
  }
}

resource "aws_launch_template" "data_refinery_compendia" {
  name = "data-refinery-compendia-${var.user}-${var.stage}"
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

  block_device_mappings {
    device_name = "/dev/sdf"
    ebs {
      volume_type = "st1"
      volume_size = 1000
      encrypted = true
      delete_on_termination = true
    }
  }

  user_data = base64encode(var.data_refinery_worker_user_data)

  tag_specifications {
      resource_type = "instance"
      tags = var.default_tags
  }

  tag_specifications {
    resource_type = "volume"
    tags = var.default_tags
  }
}
