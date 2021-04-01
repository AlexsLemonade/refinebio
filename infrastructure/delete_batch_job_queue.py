import os
from time import sleep

import boto3

AWS_REGION = os.environ["AWS_REGION"]
AWS_BATCH_QUEUE_NAME = os.environ["REFINEBIO_JOB_QUEUE_NAME"]


batch = boto3.client("batch", region_name=AWS_REGION)

batch.update_job_queue(jobQueue=AWS_BATCH_QUEUE_NAME, state="DISABLED")

while True:
    job_queues = batch.describe_job_queues(jobQueues=[AWS_BATCH_QUEUE_NAME])
    if "jobQueues" in job_queues:
        job_queue = job_queues["jobQueues"][0]
        if job_queue["state"] == "DISABLED" and job_queue["status"] != "UPDATING":
            break
    else:
        print("Unexpected response while deleting job queue.")
        print(job_queues)
        break

    sleep(3)

batch.delete_job_queue(jobQueue=AWS_BATCH_QUEUE_NAME)
