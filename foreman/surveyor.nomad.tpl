job "SURVEYOR" {
  datacenters = ["dc1"]

  type = "batch"

  parameterized {
    payload       = "forbidden"
    meta_required = [ "COMMAND", "FILE"]
  }

  group "jobs" {
    restart {
      attempts = 0
      mode = "fail"
      # delay    = "30s"
    }

    task "surveyor" {
      driver = "docker"

      # This env will be passed into the container for the job.
      env {
        AWS_ACCESS_KEY_ID = "${{AWS_ACCESS_KEY_ID}}"
        AWS_SECRET_ACCESS_KEY = "${{AWS_SECRET_ACCESS_KEY}}"
        DJANGO_SECRET_KEY = "${{DJANGO_SECRET_KEY}}"
        DJANGO_DEBUG = "${{DJANGO_DEBUG}}"

        DATABASE_NAME = "${{DATABASE_NAME}}"
        DATABASE_USER = "${{DATABASE_USER}}"
        DATABASE_PASSWORD = "${{DATABASE_PASSWORD}}"
        DATABASE_HOST = "${{DATABASE_HOST}}"
        DATABASE_PORT = "${{DATABASE_PORT}}"
        DATABASE_TIMEOUT = "${{DATABASE_TIMEOUT}}"

        RUNNING_IN_CLOUD = "False"
        WATCHTOWER="False"

        USE_S3 = "${{USE_S3}}"
        S3_BUCKET_NAME = "${{S3_BUCKET_NAME}}"
        LOCAL_ROOT_DIR = "${{LOCAL_ROOT_DIR}}"

        RAW_PREFIX = "${{RAW_PREFIX}}"
        TEMP_PREFIX = "${{RAW_PREFIX}}"
        PROCESSED_PREFIX = "${{PROCESSED_PREFIX}}"
      }

      # The resources the job will require.
      resources {
        # CPU is in AWS's CPU units.
        cpu = 500
        # Memory is in MB of RAM.
        memory = 2048
      }

      config {
        image = "${{FOREMAN_DOCKER_IMAGE}}"
        force_pull = false

        # The args to pass to the Docker container's entrypoint.
        args = [
          "${NOMAD_META_COMMAND}",
          "--file", "${NOMAD_META_FILE}",
        ]
        ${{EXTRA_HOSTS}}
        volumes = ["${{VOLUME_DIR}}:/home/user/data_store"]

        logging {
          type = "awslogs"
          config {
            awslogs-region = "${{REGION}}",
            awslogs-group = "data-refinery-log-group-${{USER}}-${{STAGE}}",
            awslogs-stream = "log-stream-nomad-docker-downloader-${{USER}}-${{STAGE}}"
          }
        }

      }
    }
  }
}
