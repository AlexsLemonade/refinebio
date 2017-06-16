from django.test import TestCase
from data_refinery_models.models import (
    SurveyJob,
    DownloaderJob,
    DownloaderJobsToBatches,
    ProcessorJob,
    ProcessorJobsToBatches,
    Batch,
    BatchStatuses
)


def get_batch():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()

    return Batch(
        survey_job=survey_job,
        source_type="ARRAY_EXPRESS",
        size_in_bytes=0,
        download_url="example.com",
        raw_format="CEL",
        processed_format="PCL",
        pipeline_required="AFFY_TO_PCL",
        platform_accession_code="A-AFFY-1",
        experiment_accession_code="E-MTAB-3050",
        experiment_title="It doesn't really matter.",
        name="CE1234.CEL",
        internal_location="A-AFFY-1/AFFY_TO_PCL/",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-05-05",
        last_uploaded_date="2017-05-05",
        status=BatchStatuses.NEW.value
    )


class DownloaderJobTestCase(TestCase):
    def test_create_job_and_relationships(self):
        """DownloaderJob, Batches, and relationships are created."""
        batches = [get_batch(), get_batch()]

        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)
        self.assertIsInstance(downloader_job.id, int)

        # If the two batch_relations got saved to the database then
        # the so must have the batches themselves so no need to test
        # that explicitly
        batch_relations = DownloaderJobsToBatches.objects.filter(
            downloader_job_id=downloader_job.id)
        self.assertEqual(len(batch_relations), 2)


class ProcessorJobTestCase(TestCase):
    def test_create_job_and_relationships(self):
        """ProcessorJob, Batches, and relationships are created."""
        batches = [get_batch(), get_batch()]

        processor_job = ProcessorJob.create_job_and_relationships(batches=batches,
                                                                  pipeline_applied="AFFY_TO_PCL")
        self.assertIsInstance(processor_job.id, int)

        # If the two batch_relations got saved to the database then
        # the so must have the batches themselves so no need to test
        # that explicitly
        batch_relations = ProcessorJobsToBatches.objects.filter(
            processor_job_id=processor_job.id)
        self.assertEqual(len(batch_relations), 2)
