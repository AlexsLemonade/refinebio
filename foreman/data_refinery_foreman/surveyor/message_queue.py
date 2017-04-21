from __future__ import absolute_import, unicode_literals
import os
from celery import Celery

"""
This module initializes a Celery app using the same name as the
Celery app defined in the workers project. This allows us to queue
tasks without needing to import the workers project.
This is desirable because the workers project has many additional
dependencies that the foreman project does not need.
"""


# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE',
                      'data_refinery_workers.settings')

app = Celery('data_refinery_workers')

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')
