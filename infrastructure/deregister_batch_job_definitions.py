import os

import boto3

AWS_REGION = os.environ["AWS_REGION"]

batch = boto3.client("batch", region_name=AWS_REGION)

# TODO: stop repeating this construction everywhere. Just set it once somewhere.
JOB_DEFINITION_PREFIX = os.environ["USER"] + "_" + os.environ["STAGE"] + "_"

# job_prefixes = [
#     "kurt-dev-",
#     "SMASHER",
# ]

# job_defs = batch.describe_job_definitions(status="ACTIVE")

# while True:
#     for job_def in job_defs["jobDefinitions"]:
#         print(f"Job: {job_def['jobDefinitionName']} {job_def['revision']}")
#         for job_prefix in job_prefixes:
#             if job_def["jobDefinitionName"].startswith(job_prefix):
#                 print(f"Deregistering job {job_def['jobDefinitionName']}")
#                 batch.deregister_job_definition(jobDefinition=job_def["jobDefinitionArn"])
#                 break

#     if "nextToken" in job_defs and job_defs["nextToken"] is not None:
#         job_defs = batch.describe_job_definitions(status="ACTIVE", nextToken=job_defs["nextToken"])
#     else:
#         break


job_definition_files = os.listdir("batch-job-templates")

job_definition_list = list(
    {JOB_DEFINITION_PREFIX + job_def.upper().split(".")[0] for job_def in job_definition_files}
)

# If we ever go over 100 definitions we'll have to make more than one call.
job_definitions = batch.describe_job_definitions(
    jobDefinitions=job_definition_list, status="ACTIVE"
)

for job_definition in job_definition_list:
    job_definitions = batch.describe_job_definitions(
        jobDefinitionName=job_definition, status="ACTIVE"
    )
    # There can be multiple revisions per job deifinition. We want them all gone.
    for job_definition_revision in job_definitions["jobDefinitions"]:
        batch.deregister_job_definition(jobDefinition=job_definition_revision["jobDefinitionArn"])
