resource "aws_security_group_rule" "data_refinery_ci_nomad" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  cidr_blocks = ["${var.host_ip}/32"]
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "aws_security_group_rule" "data_refinery_ci_postgres" {
  type = "ingress"
  from_port = 5432
  to_port = 5432
  protocol = "tcp"
  cidr_blocks = ["${var.host_ip}/32"]
  security_group_id = "${aws_security_group.data_refinery_db.id}"
}

resource "null_resource" "kill_nomad_jobs" {
  provisioner "local-exec" {
    command = <<EOF
# Wait for Nomad to get started.
start_time=$(date +%s)
diff=0
# This function checks what the status of the Nomad agent is.
# Returns an HTTP response code. i.e. 200 means everything is ok.
check_nomad_status () {
    echo $(curl --write-out %{http_code} \
                  --silent \
                  --output /dev/null \
                  http://${aws_instance.nomad_server_1.address}:4646/v1/status/leader)
}
nomad_status=$(check_nomad_status)
while [[ $diff < 300 && $nomad_status != "200" ]]; do
    sleep 1
    nomad_status=$(check_nomad_status)
    let "diff = $(date +%s) - $start_time"
done
export NOMAD_ADDR=http://${aws_instance.nomad_server_1.public_ip}:4646
# Kill Base Nomad Jobs
for job in $(nomad status | grep running | awk {'print $1'} || grep --invert-match /)
do
  nomad stop $job
done

# Kill parameterized Nomad Jobs
for job in $(nomad status | awk {'print $1'} || grep /)
do
  nomad stop $job
done
EOF
  }
  depends_on = ["aws_security_group_rule.data_refinery_ci_nomad"]
}

resource "null_resource" "run_migrations" {
  provisioner "local-exec" {
    command = <<EOF
mkdir -p migrations;
docker pull miserlou/dr_foreman:2;
docker run \
  --volume migrations \
  -e DJANGO_DEBUG=${var.django_debug} \
  -e DATABASE_NAME=${aws_db_instance.postgres_db.name} \
  -e DATABASE_HOST=${aws_db_instance.postgres_db.address} \
  -e DATABASE_USER=${var.database_user} \
  -e DATABASE_PASSWORD=${var.database_password} \
  -e DATABASE_PORT=${var.database_port} \
  -e DATABASE_TIMEOUT=${var.database_timeout} \
  -e DJANGO_SECRET_KEY=${var.django_secret_key} \
  -e RUNNING_IN_CLOUD=${var.running_in_cloud} \
  -e USE_S3=${var.use_s3} \
  -e S3_BUCKET_NAME=${var.s3_bucket_name} \
  -e LOCAL_ROOT_DIR=${var.local_root_dir} \
  -e RAW_PREFIX=${var.raw_prefix} \
  -e TEMP_PREFIX=${var.temp_prefix} \
  -e PROCESSED_PREFIX=${var.processed_prefix} \
  miserlou/dr_foreman:2 makemigrations;
docker run \
  --volume migrations \
  -e DJANGO_DEBUG=${var.django_debug} \
  -e DATABASE_NAME=${aws_db_instance.postgres_db.name} \
  -e DATABASE_HOST=${aws_db_instance.postgres_db.address} \
  -e DATABASE_USER=${var.database_user} \
  -e DATABASE_PASSWORD=${var.database_password} \
  -e DATABASE_PORT=${var.database_port} \
  -e DATABASE_TIMEOUT=${var.database_timeout} \
  -e DJANGO_SECRET_KEY=${var.django_secret_key} \
  -e RUNNING_IN_CLOUD=${var.running_in_cloud} \
  -e USE_S3=${var.use_s3} \
  -e S3_BUCKET_NAME=${var.s3_bucket_name} \
  -e LOCAL_ROOT_DIR=${var.local_root_dir} \
  -e RAW_PREFIX=${var.raw_prefix} \
  -e TEMP_PREFIX=${var.temp_prefix} \
  -e PROCESSED_PREFIX=${var.processed_prefix} \
  miserlou/dr_foreman:2 migrate auth;
docker run \
  --volume migrations \
  -e DJANGO_DEBUG=${var.django_debug} \
  -e DATABASE_NAME=${aws_db_instance.postgres_db.name} \
  -e DATABASE_HOST=${aws_db_instance.postgres_db.address} \
  -e DATABASE_USER=${var.database_user} \
  -e DATABASE_PASSWORD=${var.database_password} \
  -e DATABASE_PORT=${var.database_port} \
  -e DATABASE_TIMEOUT=${var.database_timeout} \
  -e DJANGO_SECRET_KEY=${var.django_secret_key} \
  -e RUNNING_IN_CLOUD=${var.running_in_cloud} \
  -e USE_S3=${var.use_s3} \
  -e S3_BUCKET_NAME=${var.s3_bucket_name} \
  -e LOCAL_ROOT_DIR=${var.local_root_dir} \
  -e RAW_PREFIX=${var.raw_prefix} \
  -e TEMP_PREFIX=${var.temp_prefix} \
  -e PROCESSED_PREFIX=${var.processed_prefix} \
  miserlou/dr_foreman:2 migrate;
EOF
  }
  depends_on = [
    "null_resource.kill_nomad_jobs",
    "aws_security_group_rule.data_refinery_ci_postgres"
  ]

}
