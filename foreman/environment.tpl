DJANGO_SECRET_KEY=${{DJANGO_SECRET_KEY}}
DJANGO_DEBUG=${{DJANGO_DEBUG}}
DATABASE_PORT=${{DATABASE_PORT}}
DATABASE_TIMEOUT=${{DATABASE_TIMEOUT}}
RUNNING_IN_CLOUD=${{RUNNING_IN_CLOUD}}
USE_S3=${{USE_S3}}
S3_BUCKET_NAME=${{S3_BUCKET_NAME}}
S3_COMPENDIA_BUCKET_NAME=${{S3_COMPENDIA_BUCKET_NAME}}
LOCAL_ROOT_DIR=${{LOCAL_ROOT_DIR}}
NOMAD_HOST=${{NOMAD_HOST}}
REFINEBIO_JOB_QUEUE_NAME=${{REFINEBIO_JOB_QUEUE_NAME}}
JOB_DEFINITION_PREFIX=${{USER}}_${{STAGE}}_
