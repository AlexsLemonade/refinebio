# This file contains the configuration for the database and related resources.

resource "aws_db_parameter_group" "postgres_parameters" {
  name = "postgres-parameters-${var.user}-${var.stage}"
  description = "Postgres Parameters ${var.user} ${var.stage}"
  family = "postgres11"

  parameter {
    name = "deadlock_timeout"
    value = "60000" # 60000ms = 60s
  }

  parameter {
    name = "statement_timeout"
    value = "60000" # 60000ms = 60s
  }

  parameter {
    name = "log_checkpoints"
    value = false
  }

  parameter {
    name = "work_mem"
    value = 128000
    apply_method = "pending-reboot"
  }

  # If you have a dedicated database server with 1GB or more of RAM, a
  # reasonable starting value for shared_buffers is 25% of the memory
  # in your system. There are some workloads where even larger
  # settings for shared_buffers are effective, but because PostgreSQL
  # also relies on the operating system cache, it is unlikely that an
  # allocation of more than 40% of RAM to shared_buffers will work
  # better than a smaller amount. Larger settings for shared_buffers
  # usually require a corresponding increase in max_wal_size, in order
  # to spread out the process of writing large quantities of new or
  # changed data over a longer period of time.
  parameter {
    name = "shared_buffers"
    # Note that the unit here is 8KB, so 8192 / .4 = 20480
    value = "{DBInstanceClassMemory/20480}"
    apply_method = "pending-reboot"
  }

  # https://www.2ndquadrant.com/en/blog/basics-of-tuning-checkpoints/ says:
  # Let’s also discuss the other extreme – doing very frequent
  # checkpoints (say, every second or so). That would allow keeping
  # only tiny amount of WAL and the recovery would be very fast
  # (having to replay only the tiny WAL amount). But it would also
  # turn the asynchronous writes to data files into synchronous ones,
  # seriously impacting the users (e.g. increasing COMMIT latency,
  # reducing throughput).

  # Higher WAL values means a higher recovery time, but better general
  # performance. We want to optimize for general performance over
  # recovery time.
  parameter {
    name = "max_wal_size"
    # Note that the unit here is 1MB, so this is 4GB.
    value = 4096
    apply_method = "pending-reboot"
  }

  parameter {
    name = "min_wal_size"
    # Note that the unit here is 1MB, so this is 1GB.
    value = 1024
    apply_method = "pending-reboot"
  }

  parameter {
    name = "wal_buffers"
    # Note that the unit here is 8KB so this is 16MB.
    value = 16384
    apply_method = "pending-reboot"
  }

  parameter {
    name = "effective_cache_size"
    # Note that the unit here is 8KB, so 8192 / .75 = ~10922.67
    # actually 75% of the DB instance's RAM
    value = "{DBInstanceClassMemory/10922}"
    apply_method = "pending-reboot"
  }

  parameter {
    name = "maintenance_work_mem"
    # The unit here is KB, so this is 1GB.
    value = 1048576
    apply_method = "pending-reboot"
  }

  # From https://www.postgresql.org/docs/11/runtime-config-query.html
  # Storage that has a low random read cost relative to sequential,
  # e.g., solid-state drives, might also be better modeled with a
  # lower value for random_page_cost, e.g., 1.1.
  parameter {
    name = "random_page_cost"
    value = 1.1
    apply_method = "pending-reboot"
  }

  # https://www.postgresql.org/docs/current/runtime-config-resource.html says
  # SSDs and other memory-based storage can often process many
  # concurrent requests, so the best value might be in the hundreds.
  parameter {
    name = "effective_io_concurrency"
    value = 200
    apply_method = "pending-reboot"
  }

  # These two should match how many cores we have.
  parameter {
    name = "max_worker_processes"
    value = "{DBInstanceVCPU}"
    apply_method = "pending-reboot"
  }

  parameter {
    name = "max_parallel_workers"
    value = "{DBInstanceVCPU}"
    apply_method = "pending-reboot"
  }

  parameter {
    name = "autovacuum_vacuum_scale_factor"
    value = "0.01"
    apply_method = "pending-reboot"
  }

  tags = var.default_tags
}

resource "aws_db_instance" "postgres_db" {
  identifier = "data-refinery-${var.user}-${var.stage}"
  allocated_storage = 100
  storage_type = "gp2"
  engine = "postgres"
  engine_version = "11.16"
  allow_major_version_upgrade = true
  auto_minor_version_upgrade = false
  instance_class = "db.${var.database_instance_type}"
  db_name = "data_refinery"
  port = var.database_hidden_port
  username = var.database_user
  password = var.database_password

  apply_immediately = true

  db_subnet_group_name = aws_db_subnet_group.data_refinery.name
  parameter_group_name = aws_db_parameter_group.postgres_parameters.name

  # TF is broken, but we do want this protection in prod.
  # Related: https://github.com/hashicorp/terraform/issues/5417
  # Only the prod's bucket prefix is empty.
  skip_final_snapshot = var.stage == "prod" ? false : true
  final_snapshot_identifier = var.stage == "prod" ? "data-refinery-prod-snapshot" : "none"

  enabled_cloudwatch_logs_exports = ["postgresql", "upgrade"]

  vpc_security_group_ids = [aws_security_group.data_refinery_db.id]
  multi_az = true
  publicly_accessible = true

  backup_retention_period = var.stage == "prod" ? "7" : "0"

  tags = var.default_tags
}

resource "aws_instance" "pg_bouncer" {
  ami = data.aws_ami.ubuntu.id
  instance_type = var.pg_bouncer_instance_type
  availability_zone = "${var.region}a"
  vpc_security_group_ids = [aws_security_group.data_refinery_pg.id]
  iam_instance_profile = aws_iam_instance_profile.data_refinery_api.name
  subnet_id = aws_subnet.data_refinery_1a.id
  depends_on = [aws_db_instance.postgres_db]
  key_name = aws_key_pair.data_refinery.key_name
  disable_api_termination = "false"

  # Our instance-user-data.sh script is built by Terraform at
  # apply-time so that it can put additional files onto the
  # instance. For more information see the definition of this resource.
  user_data = templatefile("workers-configuration/pg-bouncer-instance-user-data.tpl.sh",
    {
      database_host = aws_db_instance.postgres_db.address
      database_name = aws_db_instance.postgres_db.db_name
      database_password = var.database_password
      database_port = var.database_hidden_port
      database_user = var.database_user
      listen_port = var.database_port
      region = var.region
      stage = var.stage
      user = var.user
    }
  )

  tags = merge(
    var.default_tags,
    {
      Name = "pg-bouncer-${var.user}-${var.stage}"
    }
  )

  root_block_device {
    volume_type = "gp2"
    volume_size = 100

    tags = var.default_tags
  }
}
