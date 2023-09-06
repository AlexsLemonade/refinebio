import os
import time

import boto3

batch = boto3.client("batch", region_name=os.environ["AWS_REGION"])

# TODO: stop repeating this construction everywhere. Just set it once somewhere.
JOB_DEFINITION_PREFIX = os.environ["USER"] + "_" + os.environ["STAGE"] + "_"

job_names = (
    JOB_DEFINITION_PREFIX + batch_job_template.upper().split(".")[0]
    for batch_job_template in os.listdir("batch-job-templates")
)
nextToken = ""

# Have to go one by one because providing a list of job names doesn't work:
# https://github.com/boto/boto3/issues/2908
for job_name in sorted(job_names):
    while True:
        data = {
            "jobDefinitionName": job_name,
            "maxResults": 100,
            "status": "ACTIVE",
        }
        if nextToken:
            data["nextToken"] = nextToken

        response = batch.describe_job_definitions(**data)
        nextToken = response.get("nextToken", "")

        job_definitions = response.get("jobDefinitions")
        if not job_definitions:
            break

        # There can be multiple revisions per job definition. We want them all gone.
        for job_definition in job_definitions:
            batch.deregister_job_definition(jobDefinition=job_definition["jobDefinitionArn"])

        time.sleep(1)
