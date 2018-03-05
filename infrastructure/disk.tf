# Elastic File System and EFS Mount Points

resource "aws_efs_file_system" "data_refinery_efs" {
  creation_token = "data-refinery-efs-${var.user}-${var.stage}"

  tags {
    Name = "data-refinery-efs-${var.user}-${var.stage}"
  }

}

resource "aws_efs_mount_target" "data_refinery_efs_1a" {
  file_system_id = "${aws_efs_file_system.data_refinery_efs.id}"
  subnet_id = "${aws_subnet.data_refinery_1a.id}"
  security_groups = ["${aws_security_group.data_refinery_worker.id}"]
}

resource "aws_efs_mount_target" "data_refinery_efs_1b" {
  file_system_id = "${aws_efs_file_system.data_refinery_efs.id}"
  subnet_id = "${aws_subnet.data_refinery_1b.id}"
  security_groups = ["${aws_security_group.data_refinery_worker.id}"]
}
