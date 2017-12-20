job "DOWNLOADER" {
  datacenters = ["dc1"]

  type = "batch"

  parameterized {
    payload       = "forbidden"
    meta_required = ["JOB_ID"]
  }

  group "jobs" {
    restart {
      attempts = 0
      mode = "fail"
      # delay    = "30s"
    }

    task "downloader" {
      driver = "docker"

      env {
        AWS_ACCESS_KEY_ID = "REDACTED"
        AWS_SECRET_ACCESS_KEY = "REDACTED"
        DJANGO_SECRET_KEY = "THIS_IS_NOT_A_SECRET_DO_NOT_USE_IN_PROD"
        DJANGO_DEBUG = "True"

        DATABASE_NAME = "data_refinery"
        DATABASE_USER = "data_refinery_user"
        DATABASE_PASSWORD = "data_refinery_password"
        DATABASE_HOST = "database"
        DATABASE_PORT = "5432"
        DATABASE_TIMEOUT = "5"

        RUNNING_IN_CLOUD = "False"

        USE_S3 = "False"
        S3_BUCKET_NAME = "data-refinery"
        LOCAL_ROOT_DIR = "/home/user/data_store"

        RAW_PREFIX = "raw"
        TEMP_PREFIX = "temp"
        PROCESSED_PREFIX = "processed"
      }

      resources {
        cpu = 500
        memory = 2024
      }

      config {
        image = "wkurt/nomad-test:second"
        force_pull = false

        auth {
          username = "REDACTED"
          password = "REDACTED"
        }

        args = ["run_downloader_job", "--job-id", "${NOMAD_META_JOB_ID}"]

        extra_hosts = ["database:165.123.67.153"]

        volumes = ["/home/kurt/Development/data_refinery/workers/volume:/home/user/data_store"]

      }
    }
  }
}
