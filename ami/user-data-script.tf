# This is the docker apt gpg key, which we have verified is accurate
data "local_file" "docker_apt_key" {
  filename = "docker-apt-key.gpg"
}

data "template_file" "ubuntu_instance_user_data" {
  template = file("ubuntu-instance-user-data.tpl.sh")

  vars = {
    docker_apt_key = data.local_file.docker_apt_key.content
  }
}

data "template_file" "ecs_instance_user_data" {
  template = file("ecs-instance-user-data.tpl.sh")

  vars = {
    docker_apt_key = data.local_file.docker_apt_key.content
  }
}
