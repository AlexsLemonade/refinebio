job "SURVEYOR_${{RAM}}" {
  datacenters = ["dc1"]

  type = "batch"

  parameterized {
    payload       = "forbidden"
    meta_required = ["JOB_NAME", "JOB_ID"]
  }

  group "jobs" {
    restart {
      attempts = 0
      mode = "fail"
    }

    reschedule {
      attempts = 0
      unlimited = false
    }

    ephemeral_disk {
      size = "10"
    }

    task "surveyor" {
      driver = "docker"

      # This env will be passed into the container for the job.
      env {
        ${{AWS_CREDS}}
        DJANGO_SECRET_KEY = "${{DJANGO_SECRET_KEY}}"
        DJANGO_DEBUG = "${{DJANGO_DEBUG}}"

        DATABASE_NAME = "${{DATABASE_NAME}}"
        DATABASE_USER = "${{DATABASE_USER}}"
        DATABASE_PASSWORD = "${{DATABASE_PASSWORD}}"
        DATABASE_HOST = "${{DATABASE_HOST}}"
        DATABASE_PORT = "${{DATABASE_PORT}}"
        DATABASE_TIMEOUT = "${{DATABASE_TIMEOUT}}"

        RAVEN_DSN="${{RAVEN_DSN}}"
        RAVEN_DSN_API="${{RAVEN_DSN_API}}"

        RUNNING_IN_CLOUD = "False"

        USE_S3 = "${{USE_S3}}"
        S3_BUCKET_NAME = "${{S3_BUCKET_NAME}}"
        LOCAL_ROOT_DIR = "${{LOCAL_ROOT_DIR}}"
        NOMAD_HOST = "${{NOMAD_HOST}}"
        NOMAD_PORT = "${{NOMAD_PORT}}"

        LOG_LEVEL = "${{LOG_LEVEL}}"
      }

      # The resources the job will require.
      resources {
        # CPU is in AWS's CPU units.
        cpu = 500
        # Memory is in MB of RAM.
        memory = ${{RAM}}
      }

      logs {
        max_files = 1
        max_file_size = 1
      }

      config {
        image = "${{DOCKERHUB_REPO}}/${{FOREMAN_DOCKER_IMAGE}}"
        force_pull = false

        # The args to pass to the Docker container's entrypoint.
        args = [
          "python3",
          "manage.py",
          "survey_all",
          "--job-id", "${NOMAD_META_JOB_ID}",
        ]
        ${{EXTRA_HOSTS}}
        volumes = ["${{VOLUME_DIR}}:/home/user/data_store"]
        ${{LOGGING_CONFIG}}
      }
    }
  }
}
