from __future__ import absolute_import, unicode_literals
import os
from celery import Celery

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
                       related_name="salmon")
app.autodiscover_tasks(["data_refinery_workers.processors"],
                       related_name="no_op")
app.autodiscover_tasks(["data_refinery_workers.processors"],
                       related_name="transcriptome_index")
