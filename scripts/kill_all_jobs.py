import os

import boto3

AWS_REGION = os.environ["AWS_REGION"]
AWS_BATCH_QUEUE_NAME = os.environ["REFINEBIO_JOB_QUEUE_NAME"]

batch = boto3.client("batch", region_name=AWS_REGION)


for status in ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]:
    list_jobs_dict = batch.list_jobs(jobQueue=AWS_BATCH_QUEUE_NAME, jobStatus=status)

    for job in list_jobs_dict["jobSummaryList"]:
        print("Deleting job: " + job["jobId"])
        batch.terminate_job(jobId=job["jobId"], reason="kill_all_jobs")

    while "nextToken" in list_jobs_dict and list_jobs_dict["nextToken"]:
        list_jobs_dict = batch.list_jobs(
            jobQueue=AWS_BATCH_QUEUE_NAME, jobStatus=status, nextToken=list_jobs_dict["nextToken"],
        )
