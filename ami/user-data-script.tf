# This is the docker apt gpg key, which we have verified is accurate
data "local_file" "docker_apt_key" {
  filename = "docker-apt-key.gpg"
}
