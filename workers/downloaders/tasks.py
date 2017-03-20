from __future__ import absolute_import, unicode_literals
from celery import shared_task

from bioinformatics_mill_models.models import Batch, BatchStatuses, SurveyJob


@shared_task
def save_test_batch():
    print("Saving survey job")
    test_survey_job = SurveyJob(
        source_type="ARRAY_EXPRESS",
        replication_started_at="2001-04-10 15:51:24Z",
        replication_ended_at="2001-04-10 15:51:24Z",
        start_time="2001-04-10 15:51:24Z",
        end_time="2001-04-10 15:51:24Z",
    )
    test_survey_job.save()

    print("saving batch")
    test_batch = Batch(
        survey_job=test_survey_job,
        source_type="ARRAY_EXPRESS",
        size_in_bytes=1000,
        download_url="www.example.com",
        raw_format="MICRO_ARRAY",
        processed_format="PCL",
        pipeline_required="MICRO_ARRAY_TO_PCL",
        accession_code="E-MEXP-123",
        internal_location="expression_data/array_express/E-MEXP-123",
        organism=1,
        status=BatchStatuses.NEW.value
    )

    test_batch.save()

    return True
