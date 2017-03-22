from __future__ import absolute_import, unicode_literals
from celery import shared_task

from data_refinery_models.models import Batch


@shared_task
def save_test_batch():
    print("saving batch")
    test_batch = Batch(
        source_type="ARRAY_EXPRESS",
        size_in_bytes=1000,
        download_url="www.example.com",
        raw_format="MICRO_ARRAY",
        processed_format="PCL",
        processor_required=1,
        accession_code="E-MEXP-123",
        internal_location="expression_data/array_express/E-MEXP-123",
        organism=1,
        status="PROCESSED"
    )

    test_batch.save()

    return True
