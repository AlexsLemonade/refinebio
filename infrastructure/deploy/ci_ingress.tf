resource "aws_security_group_rule" "data_refinery_ci_nomad" {
  type = "ingress"
  from_port = 4646
  to_port = 4646
  protocol = "tcp"
  cidr_blocks = ["${var.HOST_IP}/32"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "aws_security_group_rule" "data_refinery_ci_postgres" {
  type = "ingress"
  from_port = 5432
  to_port = 5432
  protocol = "tcp"
  cidr_blocks = ["${var.HOST_IP}/32"]
  security_group_id = "${aws_security_group.data_refinery_worker.id}"
}

resource "null_resource" "kill_nomad_jobs" {
  provisioner "local-exec" {
    command = <<EOF
# Wait for Nomad to get started.
start_time=$(date +%s)
diff=0
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
}

resource "null_resource" "run_migrations" {
  provisioner "local-exec" {
    command = <<EOF
mkdir migrations
docker run \
  --volume migrations \
  -e DATABASE_NAME=data_refinery \
  -e DATABASE_HOST=${aws_db_instance.postgres_db.address} \
  -e DATABASE_USER=data_refinery_user \
  -e DATABASE_PASSWORD=data_refinery_password \
  -e DATABASE_PORT=${var.database_port} \
  -e DATABASE_TIMEOUT=${var.database_timeout} \
  -e DJANGO_SECRET_KEY=${var.django_secret_key} \
  miserlou/dr_foreman:latest python3.6 manage.py makemigrations && python3.6 manage.py migrate
EOF
  }
}
