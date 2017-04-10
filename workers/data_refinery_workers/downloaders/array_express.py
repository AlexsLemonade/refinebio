from __future__ import absolute_import, unicode_literals
from celery import shared_task

from data_refinery_models.models import SurveyJob


@shared_task
def download_array_express(downloader_job_id):
    test_survey_job = SurveyJob(
        source_type="JUST_TESTING",
        replication_started_at="2001-04-10 15:51:24Z",
        replication_ended_at="2001-04-10 15:51:24Z",
        start_time="2001-04-10 15:51:24Z",
        end_time="2001-04-10 15:51:24Z",
    )
    test_survey_job.save()
