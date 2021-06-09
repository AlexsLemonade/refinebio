import os
from time import sleep

import boto3

AWS_REGION = os.environ["AWS_REGION"]
AWS_BATCH_QUEUE_ALL_NAMES = os.environ["REFINEBIO_JOB_QUEUE_ALL_NAMES"].split(",")

batch = boto3.client("batch", region_name=AWS_REGION)

# First disable each job queue.
for batch_queue_name in AWS_BATCH_QUEUE_ALL_NAMES:
    try:
        batch.update_job_queue(jobQueue=batch_queue_name, state="DISABLED")
    except Exception as e:
        # If the job queue doesn't exist, that's cool, we were trying to delete it anyway.
        pass

# Then wait for each one to be disabled so it can be deleted.
for batch_queue_name in AWS_BATCH_QUEUE_ALL_NAMES:
    while True:
        job_queues = batch.describe_job_queues(jobQueues=[batch_queue_name])
        if "jobQueues" in job_queues:
            job_queue = job_queues["jobQueues"][0]
            if job_queue["state"] == "DISABLED" and job_queue["status"] != "UPDATING":
                break
        else:
            print(f"Unexpected response while describing job queue {batch_queue_name}.")
            break

        sleep(3)

    batch.delete_job_queue(jobQueue=batch_queue_name)
