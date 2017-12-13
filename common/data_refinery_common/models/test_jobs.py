from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    ProcessorJob,
    Batch,
    BatchStatuses,
    File
)


def get_batch():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="ARRAY_EXPRESS",
        pipeline_required="AFFY_TO_PCL",
        platform_accession_code="A-AFFY-1",
        experiment_accession_code="E-MTAB-3050",
        experiment_title="It doesn't really matter.",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-05-05",
        last_uploaded_date="2017-05-05",
        status=BatchStatuses.NEW.value
    )
    batch.save()

    File(
        size_in_bytes=0,
        download_url="example.com",
        raw_format="CEL",
        processed_format="PCL",
        name="CE1234.CEL",
        internal_location="A-AFFY-1/AFFY_TO_PCL/",
        batch=batch
    ).save()
    return batch


class DownloaderJobTestCase(TestCase):
    def test_create_job_and_relationships(self):
        """DownloaderJob, Batches, and relationships are created."""
        batches = [get_batch(), get_batch()]

        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches,
                                                                    downloader_task="test")
        self.assertIsInstance(downloader_job.id, int)

        batches_for_job = downloader_job.batches.all()
        self.assertEqual(len(batches_for_job), 2)

        self.assertEqual(batches[0].downloaderjob_set.get(), downloader_job)


class ProcessorJobTestCase(TestCase):
    def test_create_job_and_relationships(self):
        """ProcessorJob, Batches, and relationships are created."""
        batches = [get_batch(), get_batch()]

        processor_job = ProcessorJob.create_job_and_relationships(batches=batches,
                                                                  pipeline_applied="AFFY_TO_PCL")
        self.assertIsInstance(processor_job.id, int)

        batches_for_job = processor_job.batches.all()
        self.assertEqual(len(batches_for_job), 2)

        self.assertEqual(batches[0].processorjob_set.get(), processor_job)
