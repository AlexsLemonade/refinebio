# This script installs Nomad.
data "local_file" "install_nomad_script" {
  filename = "../scripts/install_nomad.sh"
}

# This is the docker apt gpg key, which we have verified is accurate
data "local_file" "docker_apt_key" {
  filename = "docker-apt-key.gpg"
}

data "template_file" "instance_user_data" {
  template = file("instance-user-data.tpl.sh")

  vars {
    install_nomad_script = data.local_file.install_nomad_script.content
    docker_apt_key = data.local_file.docker_apt_key.content
  }
}
