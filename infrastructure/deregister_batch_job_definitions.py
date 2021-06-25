import os

import boto3

AWS_REGION = os.environ["AWS_REGION"]

batch = boto3.client("batch", region_name=AWS_REGION)

# TODO: stop repeating this construction everywhere. Just set it once somewhere.
JOB_DEFINITION_PREFIX = os.environ["USER"] + "_" + os.environ["STAGE"] + "_"

job_definition_files = os.listdir("batch-job-templates")

job_definition_list = list(
    {JOB_DEFINITION_PREFIX + job_def.upper().split(".")[0] for job_def in job_definition_files}
)

# Have to go one by one because providing a list of job names doesn't work:
# https://github.com/boto/boto3/issues/2908
for job_definition in job_definition_list:
    job_definitions = batch.describe_job_definitions(
        jobDefinitionName=job_definition, status="ACTIVE"
    )
    # There can be multiple revisions per job deifinition. We want them all gone.
    for job_definition_revision in job_definitions["jobDefinitions"]:
        batch.deregister_job_definition(jobDefinition=job_definition_revision["jobDefinitionArn"])
