from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
import nomad

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE',
                      'data_refinery_workers.settings')

app = Celery('data_refinery_workers')

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

app.autodiscover_tasks(["data_refinery_workers.downloaders"],
                       related_name="array_express")
app.autodiscover_tasks(["data_refinery_workers.downloaders"],
                       related_name="sra")
app.autodiscover_tasks(["data_refinery_workers.downloaders"],
                       related_name="transcriptome_index")
app.autodiscover_tasks(["data_refinery_workers.processors"],
                       related_name="array_express")
app.autodiscover_tasks(["data_refinery_workers.processors"],
                       related_name="no_op")
app.autodiscover_tasks(["data_refinery_workers.processors"],
                       related_name="transcriptome_index")


def send_job(job_name: str, job_id: int):
    nomad_client = nomad.Nomad("database", timeout=5)
    nomad_client.job.dispatch_job(job_name, meta={"JOB_ID": str(job_id)})
