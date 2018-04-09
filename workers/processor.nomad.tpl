# This job name is temporary. Once resource requirements for various
# jobs have been properly estimated there should be different job
# specifications for each processor type.
job "PROCESSOR" {
  datacenters = ["dc1"]

  type = "batch"

  parameterized {
    payload       = "forbidden"
    meta_required = [ "JOB_NAME", "JOB_ID"]
  }

  group "jobs" {
    restart {
      attempts = 0
      mode = "fail"
      # delay    = "30s"
    }

    task "processor" {
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

        RUNNING_IN_CLOUD = "${{RUNNING_IN_CLOUD}}"

        USE_S3 = "${{USE_S3}}"
        S3_BUCKET_NAME = "${{S3_BUCKET_NAME}}"
        LOCAL_ROOT_DIR = "${{LOCAL_ROOT_DIR}}"

        RAW_PREFIX = "${{RAW_PREFIX}}"
        TEMP_PREFIX = "${{RAW_PREFIX}}"
        PROCESSED_PREFIX = "${{PROCESSED_PREFIX}}"

        NOMAD_HOST = "${{NOMAD_HOST}}"
      }

      # The resources the job will require.
      resources {
        # CPU is in AWS's CPU units.
        cpu = 1024
        # Memory is in MB of RAM.
        memory = 4048
      }

      config {
        image = "${{WORKERS_DOCKER_IMAGE}}"
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
