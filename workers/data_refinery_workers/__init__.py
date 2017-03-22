from __future__ import absolute_import, unicode_literals

# This will make sure the app is always imported when
# Django starts so that shared_task will use this app.
from .task_runner import app as task_runner

__all__ = ['task_runner']
