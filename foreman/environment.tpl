AWS_ACCESS_KEY_ID=${{AWS_ACCESS_KEY_ID_WORKER}}
AWS_SECRET_ACCESS_KEY=${{AWS_SECRET_ACCESS_KEY_WORKER}}
DJANGO_SECRET_KEY=${{DJANGO_SECRET_KEY}}
DJANGO_DEBUG=${{DJANGO_DEBUG}}
DATABASE_PORT=${{DATABASE_PORT}}
DATABASE_TIMEOUT=15
RUNNING_IN_CLOUD=True
USE_S3=${{USE_S3}}
S3_BUCKET_NAME=${{S3_BUCKET_NAME}}
LOCAL_ROOT_DIR=${{LOCAL_ROOT_DIR}}
NOMAD_HOST=${{NOMAD_HOST}}
RAW_PREFIX=raw
TEMP_PREFIX=temp
PROCESSED_PREFIX=processed

