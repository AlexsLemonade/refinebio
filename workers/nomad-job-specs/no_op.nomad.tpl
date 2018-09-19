job "NO_OP_${{INDEX}}_${{RAM}}" {
  datacenters = ["dc1"]

  type = "batch"
  priority = 50

  parameterized {
    payload       = "forbidden"
    meta_required = ["JOB_NAME", "JOB_ID"]
  }

  group "jobs" {
    restart {
      attempts = 0
      mode = "fail"
      # delay    = "30s"
    }

    ephemeral_disk {
      size    = "1"
    }

    task "no_op" {
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

        NOMAD_HOST = "${{NOMAD_HOST}}"
        NOMAD_PORT = "${{NOMAD_PORT}}"

        RUNNING_IN_CLOUD = "${{RUNNING_IN_CLOUD}}"

        USE_S3 = "${{USE_S3}}"
        S3_BUCKET_NAME = "${{S3_BUCKET_NAME}}"
        LOCAL_ROOT_DIR = "${{LOCAL_ROOT_DIR}}"
        EBS_INDEX = "${{INDEX}}"

      }

      # The resources the job will require.
      resources {
        # CPU is in AWS's CPU units.
        cpu = 256
        # Memory is in MB of RAM.
        memory = ${{RAM}}
      }

      logs {
        max_files = 1
        max_file_size = 1
      }

      constraint {
        attribute = "${meta.volume_index}"
        operator  = "="
        value     = "${{INDEX}}"
      }

      config {
        image = "${{DOCKERHUB_REPO}}/${{NO_OP_DOCKER_IMAGE}}"
        force_pull = false

        # The args to pass to the Docker container's entrypoint.
        args = [
          "python3",
          "manage.py",
          "run_processor_job",
          "--job-name", "${NOMAD_META_JOB_NAME}",
          "--job-id", "${NOMAD_META_JOB_ID}"
        ]
        ${{EXTRA_HOSTS}}
        volumes = ["${{VOLUME_DIR}}:/home/user/data_store"]
        ${{LOGGING_CONFIG}}
      }
    }
  }
}
