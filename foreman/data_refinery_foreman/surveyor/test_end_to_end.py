import glob
import logging
import os
import shutil
import time
from datetime import datetime, timedelta
from test.support import EnvironmentVarGuard  # Python >=3
from unittest import skip
from unittest.mock import patch

from django.test import TransactionTestCase, tag
from django.utils import timezone

import pandas as pd
import scipy.stats
import vcr

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_foreman.foreman.main import retry_lost_downloader_jobs
from data_refinery_foreman.surveyor import surveyor
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logging.getLogger("vcr").setLevel(logging.WARN)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
CASSETTES_DIR = "/home/user/data_store/cassettes/"
LOOP_TIME = 5  # seconds
MAX_WAIT_TIME = timedelta(minutes=60)

# We don't want our end-to-end tests to be dependent upon external
# services, so we host the files we would download fromt hem on S3. We
# then have to replace the source_urls of downloader jobs that are
# generated for end-to-end tests to pull from S3 instead of the
# original location. The following constants and function wait until
# the original files have been created and then swap out their
# source_urls based on the mapping below.
# We have to save references to the actual surveyor classes before
# they get overwritten with mocks.
ORIGINAL_ARRAY_EXPRESS_SURVEYOR = surveyor.ArrayExpressSurveyor
ORIGINAL_SRA_SURVEYOR = surveyor.SraSurveyor
ORIGINAL_TRANSCRIPTOME_SURVEYOR = surveyor.TranscriptomeIndexSurveyor
ORIGINAL_GEO_SURVEYOR = surveyor.GeoSurveyor

EXTERNAL_FILE_URL_MAPPING = {
    # Transcriptome:
    "ftp://ftp.ensembl.org/pub/release-99/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.99.gtf.gz": "https://data-refinery-test-assets.s3.amazonaws.com/end_to_end_downloads/Caenorhabditis_elegans.WBcel235.99.gtf.gz",  # noqa
    "ftp://ftp.ensembl.org/pub/release-99/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz": "https://data-refinery-test-assets.s3.amazonaws.com/end_to_end_downloads/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz",  # noqa
    # No Op:
    "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-3303/E-GEOD-3303.processed.1.zip": "https://data-refinery-test-assets.s3.amazonaws.com/end_to_end_downloads/E-GEOD-3303.processed.1.zip",  # noqa
    # GEO:
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102571/miniml/GSE102571_family.xml.tgz": "https://data-refinery-test-assets.s3.amazonaws.com/end_to_end_downloads/GSE102571_family.xml.tgz",  # noqa
}


def build_surveyor_init_mock(source_type):
    if source_type == "ARRAY_EXPRESS":
        original_surveyor = ORIGINAL_ARRAY_EXPRESS_SURVEYOR
    if source_type == "SRA":
        original_surveyor = ORIGINAL_SRA_SURVEYOR
    if source_type == "TRANSCRIPTOME_INDEX":
        original_surveyor = ORIGINAL_TRANSCRIPTOME_SURVEYOR
    if source_type == "GEO":
        original_surveyor = ORIGINAL_GEO_SURVEYOR

    def mock_init_surveyor(survey_job):
        ret_value = original_surveyor(survey_job)

        def mock_queue_downloader_job_for_original_files(
            original_files, experiment_accession_code: str = None, is_transcriptome: bool = False,
        ):
            for original_file in original_files:
                original_file.source_url = EXTERNAL_FILE_URL_MAPPING[original_file.source_url]
                original_file.save()

            return original_surveyor.queue_downloader_job_for_original_files(
                ret_value, original_files, experiment_accession_code, is_transcriptome
            )

        def mock_queue_downloader_jobs(experiment, samples):
            # We don't want to change the same file's url more than
            # once and sometimes multiple samples are associated with
            # the same file.
            original_file_ids = set()
            for sample in samples:
                for original_file in sample.original_files.all():
                    if original_file.id not in original_file_ids:
                        try:
                            original_file.source_url = EXTERNAL_FILE_URL_MAPPING[
                                original_file.source_url
                            ]
                        except KeyError as e:
                            log_message = (
                                "The tests attempted to access a URL that is not in"
                                " EXTERNAL_FILE_URL_MAPPING. This is most likely because you've"
                                " added a test that mocks the surveyor so that the DownloaderJobs"
                                " download from S3 instead of the external service. To fix this"
                                " you should download the file, upload it to S3, and then add a"
                                " mapping from its original URL to its URL in S3."
                            )
                            logger.warn(log_message)
                            raise e

                        original_file.save()
                        original_file_ids.add(original_file.id)

            return original_surveyor.queue_downloader_jobs(ret_value, experiment, samples)

        ret_value.queue_downloader_job_for_original_files = (
            mock_queue_downloader_job_for_original_files
        )
        ret_value.queue_downloader_jobs = mock_queue_downloader_jobs
        return ret_value

    return mock_init_surveyor


def wait_for_job(job, job_class: type, start_time: datetime, loop_time: int = None):
    """Monitors the `job_class` table for when `job` is done."""
    loop_time = loop_time if loop_time else LOOP_TIME
    last_log_time = timezone.now()
    while job.success is None and timezone.now() - start_time < MAX_WAIT_TIME:
        time.sleep(loop_time)

        # Don't log statuses more often than every 15 seconds.
        if timezone.now() - last_log_time > timedelta(seconds=15):
            logger.info("Still polling the %s with ID %s.", job_class.__name__, job.nomad_job_id)
            last_log_time = timezone.now()

        job.refresh_from_db()

    if timezone.now() - start_time > MAX_WAIT_TIME:
        logger.error("%s job timed out!", job_class.__name__)

    return job


# TransactionTestCase makes database calls complete before the test
# ends.  Otherwise the workers wouldn't actually be able to find the
# job in the database because it'd be stuck in a transaction.
#
# Unfortunately, it's more unreliable and tends to leave things in the
# database, so let's manually clear it before every test
class EndToEndTestCase(TransactionTestCase):
    def setUp(self):
        Experiment.objects.all().delete()
        ExperimentAnnotation.objects.all().delete()
        ExperimentSampleAssociation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()
        OriginalFile.objects.all().delete()
        OriginalFileSampleAssociation.objects.all().delete()
        SampleResultAssociation.objects.all().delete()
        ComputationalResult.objects.all().delete()
        ComputationalResultAnnotation.objects.all().delete()
        SampleComputedFileAssociation.objects.all().delete()
        ComputedFile.objects.all().delete()
        DownloaderJob.objects.all().delete()
        DownloaderJobOriginalFileAssociation.objects.all().delete()
        ProcessorJob.objects.all().delete()
        ProcessorJobOriginalFileAssociation.objects.all().delete()


class NoOpEndToEndTestCase(EndToEndTestCase):
    @tag("slow")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.no_op.yaml"), ignore_hosts=["nomad"],
    )
    @patch("data_refinery_foreman.surveyor.surveyor.ArrayExpressSurveyor")
    def test_no_op(self, mock_surveyor):
        """Survey, download, then process an experiment we know is NO_OP."""

        mock_surveyor.side_effect = build_surveyor_init_mock("ARRAY_EXPRESS")

        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            accession_code = "E-GEOD-3303"
            survey_job = surveyor.survey_experiment(accession_code, "ARRAY_EXPRESS")

            self.assertTrue(survey_job.success)

            downloader_jobs = DownloaderJob.objects.all()
            self.assertGreater(downloader_jobs.count(), 0)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            for downloader_job in downloader_jobs:
                downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time)
                self.assertTrue(downloader_job.success)

            processor_jobs = ProcessorJob.objects.all().exclude(
                abort=True
            )  # exclude aborted processor jobs
            self.assertGreater(processor_jobs.count(), 0)

            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            start_time = timezone.now()
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                if not processor_job.success:
                    logger.error(processor_job.failure_reason)
                self.assertTrue(processor_job.success)

            # Test that the unsurveyor deletes all objects related to the experiment
            purge_experiment(accession_code)

            self.assertEqual(Experiment.objects.all().count(), 0)
            self.assertEqual(ExperimentAnnotation.objects.all().count(), 0)
            self.assertEqual(ExperimentSampleAssociation.objects.all().count(), 0)
            self.assertEqual(Sample.objects.all().count(), 0)
            self.assertEqual(SampleAnnotation.objects.all().count(), 0)
            self.assertEqual(OriginalFile.objects.all().count(), 0)
            self.assertEqual(OriginalFileSampleAssociation.objects.all().count(), 0)
            self.assertEqual(SampleResultAssociation.objects.all().count(), 0)
            self.assertEqual(ComputationalResult.objects.all().count(), 0)
            self.assertEqual(ComputationalResultAnnotation.objects.all().count(), 0)
            self.assertEqual(SampleComputedFileAssociation.objects.all().count(), 0)
            self.assertEqual(ComputedFile.objects.all().count(), 0)
            self.assertEqual(DownloaderJob.objects.all().count(), 0)
            self.assertEqual(DownloaderJobOriginalFileAssociation.objects.all().count(), 0)
            self.assertEqual(ProcessorJob.objects.all().count(), 0)
            self.assertEqual(ProcessorJobOriginalFileAssociation.objects.all().count(), 0)


class ArrayexpressRedownloadingTestCase(EndToEndTestCase):
    @tag("slow")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.array_express_redownloading.yaml"),
        ignore_hosts=["nomad"],
    )
    @patch("data_refinery_foreman.surveyor.surveyor.ArrayExpressSurveyor")
    def test_array_express_redownloading(self, mock_surveyor):
        """Survey, download, then process an experiment we know is NO_OP."""

        mock_surveyor.side_effect = build_surveyor_init_mock("ARRAY_EXPRESS")
        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            NUM_SAMPLES_IN_EXPERIMENT = 12
            accession_code = "E-GEOD-3303"
            survey_job = surveyor.survey_experiment(accession_code, "ARRAY_EXPRESS")

            self.assertTrue(survey_job.success)

            # All of this experiment's samples are contained in the
            # same archive, so only one job is needed.
            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 1)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            # We want to try and delete the file as quickly as
            # possible, so pass a short loop time and let the waiting
            # loop spin really fast so we lose as little time as
            # possible.
            downloader_job = wait_for_job(downloader_jobs[0], DownloaderJob, start_time, 0.1)
            self.assertTrue(downloader_job.success)

            # Now we're going to delete one of the extracted files but not the other.
            deleted_file = OriginalFile.objects.filter(is_archive=False).first()
            self.assertIsNotNone(deleted_file)
            deleted_file.delete_local_file()

            # The one downloader job should have extracted all the files
            # and created as many processor jobs.
            processor_jobs = ProcessorJob.objects.all()
            self.assertEqual(processor_jobs.count(), NUM_SAMPLES_IN_EXPERIMENT)

            doomed_processor_job = deleted_file.processor_jobs.all()[0]
            logger.info(
                "Waiting on processor Nomad job %s to fail because it realized it is missing a file.",
                doomed_processor_job.nomad_job_id,
            )

            start_time = timezone.now()
            doomed_processor_job = wait_for_job(doomed_processor_job, ProcessorJob, start_time)
            self.assertTrue(doomed_processor_job.abort)

            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should now be two.
            downloader_jobs = DownloaderJob.objects.all().order_by("-id")
            self.assertEqual(downloader_jobs.count(), 2)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info("Waiting on downloader Nomad job %s", recreated_job.nomad_job_id)
            recreated_job = wait_for_job(recreated_job, DownloaderJob, start_time)
            self.assertTrue(recreated_job.success)

            # Once the Downloader job succeeds, it should create one
            # and only one processor job, after which the total goes back up
            # to NUM_SAMPLES_IN_EXPERIMENT:
            processor_jobs = ProcessorJob.objects.all().exclude(
                abort=True
            )  # exclude aborted processor jobs
            logger.error(processor_jobs)
            self.assertEqual(processor_jobs.count(), NUM_SAMPLES_IN_EXPERIMENT)

            # And finally we can make sure that all of the
            # processor jobs were successful, including the one that
            # got recreated.
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                self.assertTrue(processor_job.success)


class GeoArchiveRedownloadingTestCase(EndToEndTestCase):
    @tag("slow")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.geo_archive_redownloading.yaml"),
        ignore_hosts=["nomad"],
    )
    def test_geo_archive_redownloading(self):
        """Survey, download, then process an experiment we know is NO_OP.

        All the data for the experiment are in the same archive, which
        is one of ways we expect GEO data to come.

        This is another test which uses Aspera so it unfortunately
        cannot be made to run without relying on NCBI's aspera server.
        """
        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            accession_code = "GSE102571"
            survey_job = surveyor.survey_experiment(accession_code, "GEO")

            self.assertTrue(survey_job.success)

            # This experiment has multiple samples that are contained in the
            # same archive, so only one job is needed.
            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 1)

            logger.info(
                "Survey Job finished, waiting for Downloader Job with Nomad ID %s to complete.",
                downloader_jobs[0].nomad_job_id,
            )

            # We're going to spin as fast as we can so we can delete
            # the file in between when the downloader job finishes and
            # the processor job starts.
            start_time = timezone.now()
            file_deleted = False
            while not file_deleted and timezone.now() - start_time < MAX_WAIT_TIME:
                non_archive_files = OriginalFile.objects.filter(is_archive=False)
                for original_file in non_archive_files:
                    if original_file.absolute_file_path and os.path.exists(
                        original_file.absolute_file_path
                    ):
                        os.remove(original_file.absolute_file_path)
                        file_deleted = True
                        break

            downloader_job = wait_for_job(downloader_jobs[0], DownloaderJob, start_time)
            self.assertTrue(downloader_job.success)

            try:
                doomed_processor_job = original_file.processor_jobs.all()[0]
            except Exception:
                # The doomed job may be aborted before we can get
                # it. This is fine, we just can't look at it.
                doomed_processor_job = None

            if doomed_processor_job:
                logger.info(
                    "Waiting on processor Nomad job %s to fail because it realized it is missing a file.",
                    doomed_processor_job.nomad_job_id,
                )

                start_time = timezone.now()
                doomed_processor_job = wait_for_job(doomed_processor_job, ProcessorJob, start_time)
                self.assertTrue(doomed_processor_job.abort)

            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should now be two.
            downloader_jobs = DownloaderJob.objects.all().order_by("-id")
            self.assertEqual(downloader_jobs.count(), 2)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info("Waiting on downloader Nomad job %s", recreated_job.nomad_job_id)
            recreated_job = wait_for_job(recreated_job, DownloaderJob, start_time)
            self.assertTrue(recreated_job.success)

            # And finally we can make sure that all of the
            # processor jobs were successful, including the one that
            # got recreated.
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            processor_jobs = ProcessorJob.objects.all().exclude(
                abort=True
            )  # exclude aborted processor jobs
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                if not processor_job.success:
                    logger.error(processor_job.failure_reason)
                self.assertTrue(processor_job.success)

            # Apparently this experiment has a variable number of
            # files because GEO processed experiments sometimes do...
            # However this is okay because there's at least one file
            # per sample, so each sample will get processed at least
            # once and it's the best we can do with the state of GEO.
            # Anyway, all of that is an explanation for why we count
            # how many samples there are rather than just expecting
            # how many we know the experiment has.
            self.assertEqual(processor_jobs.count(), Sample.objects.all().count())


class GeoCelgzRedownloadingTestCase(EndToEndTestCase):
    @tag("slow")
    @tag("affymetrix")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.geo_celgz_redownloading.yaml"),
        ignore_hosts=["nomad"],
    )
    def test_geo_celgz_redownloading(self):
        """Survey, download, then process an experiment we know is Affymetrix.

        Each of the experiment's samples are in their own .cel.gz
        file, which is another way we expect GEO data to come.

        This is another test which uses Aspera so it unfortunately
        cannot be made to run without relying on NCBI's aspera server.
        """
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            # Clear out pre-existing work dirs so there's no conflicts:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="DANIO_RERIO", taxonomy_id=7955, is_scientific_name=True)
            organism.save()

            accession_code = "GSE8724"
            survey_job = surveyor.survey_experiment(accession_code, "GEO")

            SAMPLES_IN_EXPERIMENT = 3

            self.assertTrue(survey_job.success)

            # This experiment's samples each have their own file so
            # they each get their own downloader job.
            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), SAMPLES_IN_EXPERIMENT)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")

            # We're going to spin as fast as we can so we can delete
            # the file in between when the downloader jobs finishes and
            # the processor job starts.
            start_time = timezone.now()
            file_deleted = False
            while not file_deleted and timezone.now() - start_time < MAX_WAIT_TIME:
                non_archive_files = OriginalFile.objects.filter(is_archive=False)
                for original_file in non_archive_files:
                    if original_file.absolute_file_path and os.path.exists(
                        original_file.absolute_file_path
                    ):
                        os.remove(original_file.absolute_file_path)
                        file_deleted = True
                        break

            # Wait for each of the DownloaderJobs to finish
            for downloader_job in downloader_jobs:
                downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time)
                self.assertTrue(downloader_job.success)

            try:
                doomed_processor_job = original_file.processor_jobs.all()[0]
            except Exception:
                # The doomed job may be aborted before we can get
                # it. This is fine, we just can't look at it.
                doomed_processor_job = None

            if doomed_processor_job:
                logger.info(
                    "Waiting on processor Nomad job %s to fail because it realized it is missing a file.",
                    doomed_processor_job.nomad_job_id,
                )

                start_time = timezone.now()
                doomed_processor_job = wait_for_job(doomed_processor_job, ProcessorJob, start_time)
                self.assertTrue(doomed_processor_job.abort)

            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should
            # now be SAMPLES_IN_EXPERIMENT + 1 downloader jobs.
            downloader_jobs = DownloaderJob.objects.all().order_by("-id")
            self.assertEqual(downloader_jobs.count(), SAMPLES_IN_EXPERIMENT + 1)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info("Waiting on downloader Nomad job %s", recreated_job.nomad_job_id)
            recreated_job = wait_for_job(recreated_job, DownloaderJob, start_time)
            self.assertTrue(recreated_job.success)

            # And finally we can make sure that all of the processor
            # jobs were successful, including the one that got
            # recreated. The processor job that recreated that job has
            # abort=True
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            processor_jobs = ProcessorJob.objects.all().exclude(abort=True)  # exclude aborted jobs
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                self.assertTrue(processor_job.success)

            self.assertEqual(processor_jobs.count(), SAMPLES_IN_EXPERIMENT)


class TranscriptomeRedownloadingTestCase(EndToEndTestCase):
    @tag("slow")
    @tag("transcriptome")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.transcriptome_redownloading.yaml"),
        ignore_hosts=["nomad"],
    )
    @patch("data_refinery_foreman.surveyor.surveyor.TranscriptomeIndexSurveyor")
    def test_transcriptome_redownloading(self, mock_surveyor):
        """Survey, download, then process a transcriptome index."""

        mock_surveyor.side_effect = build_surveyor_init_mock("TRANSCRIPTOME_INDEX")

        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            # I'm not sure why, but sometimes there are already downloader jobs
            # in the database from previous tests even though they should be
            # removed, so pause a bit
            time.sleep(10)
            downloader_jobs = DownloaderJob.objects.all()
            for job in downloader_jobs:
                print(job)
                print(job.accession_code)
            self.assertEqual(downloader_jobs.count(), 0)

            for length in ["LONG", "SHORT"]:
                work_dir_glob = (
                    LOCAL_ROOT_DIR + "/Caenorhabditis_elegans/" + length + "/processor_job_*"
                )
                for work_dir in glob.glob(work_dir_glob):
                    shutil.rmtree(work_dir)

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(
                name="CAENORHABDITIS_ELEGANS", taxonomy_id=6239, is_scientific_name=True
            )
            organism.save()

            # Make sure that we can delete the file before the processors begin
            # by preventing the downloaders from sending the processors
            # automatically. We send the jobs manually later
            no_dispatch = EnvironmentVarGuard()
            no_dispatch.set("AUTO_DISPATCH_NOMAD_JOBS", "False")
            with no_dispatch:
                survey_job = surveyor.survey_transcriptome_index(
                    "Caenorhabditis elegans", "Ensembl"
                )

            self.assertTrue(survey_job.success)

            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 1)

            logger.info(
                "Survey Job finished, waiting for Downloader Job with Nomad ID %s to complete.",
                downloader_jobs[0].nomad_job_id,
            )

            downloader_job = wait_for_job(downloader_jobs[0], DownloaderJob, timezone.now())
            self.assertTrue(downloader_job.success)

            og_file_to_delete = OriginalFile.objects.all()[0]
            os.remove(og_file_to_delete.absolute_file_path)

            processor_jobs = ProcessorJob.objects.all()
            for processor_job in processor_jobs:
                # FIXME: we run these in serial because of
                # https://github.com/AlexsLemonade/refinebio/issues/2321
                send_job(
                    ProcessorPipeline[processor_job.pipeline_applied],
                    job=processor_job,
                    is_dispatch=True,
                )
                try:
                    wait_for_job(processor_job, ProcessorJob, timezone.now())
                except Exception:
                    pass

            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should now be two.
            downloader_jobs = DownloaderJob.objects.all().order_by("-id")
            self.assertEqual(downloader_jobs.count(), 2)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info("Waiting on downloader Nomad job %s", recreated_job.nomad_job_id)
            recreated_job = wait_for_job(recreated_job, DownloaderJob, timezone.now())
            self.assertTrue(recreated_job.success)

            # Once the Downloader job succeeds, it should create two
            # processor jobs, one for long and one for short indices.:
            processor_jobs = ProcessorJob.objects.all()
            self.assertEqual(processor_jobs.count(), 4)

            # Wait for the processor jobs to be dispatched
            time.sleep(15)

            # And finally we can make sure that both of the
            # processor jobs were successful, including the one that
            # got recreated.
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            successful_processor_jobs = []
            for processor_job in processor_jobs:
                processor_job.refresh_from_db()
                # One of the calls to wait_for_job will fail if the
                # job aborts before it we selected all the
                # processor jobs.
                processor_job = wait_for_job(processor_job, ProcessorJob, timezone.now())
                if processor_job.success:
                    successful_processor_jobs.append(processor_job)

            # While one of the original ProcessorJobs will  be aborted
            # it is hard to be sure of what will happen
            # to the other because of the racing that happens between
            # processor jobs getting started and us deleting the files
            # they need.
            # Therefore, we're just going to verify that one processor
            # job completed successfully for each length, since that
            # is the main thing we need.
            has_long = False
            has_short = False
            for processor_job in successful_processor_jobs:
                if processor_job.pipeline_applied == "TRANSCRIPTOME_INDEX_LONG":
                    has_long = True
                elif processor_job.pipeline_applied == "TRANSCRIPTOME_INDEX_SHORT":
                    has_short = True

            self.assertTrue(has_long)
            self.assertTrue(has_short)

            pj = ProcessorJob()
            pj.pipeline_applied = "SALMON"
            pj.ram_amount = 4096
            pj.save()

            c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS", taxonomy_id=9606)

            samp = Sample()
            samp.accession_code = "SALMON"  # So the test files go to the right place
            samp.organism = c_elegans
            samp.source_database = "SRA"
            samp.technology = "RNA-SEQ"
            samp.save()

            og_file = OriginalFile()
            filename = "ERR1562482.sra"
            og_file.source_filename = filename
            og_file.filename = filename
            og_file.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/" + filename
            og_file.is_downloaded = True
            og_file.save()

            og_file_samp_assoc = OriginalFileSampleAssociation()
            og_file_samp_assoc.original_file = og_file
            og_file_samp_assoc.sample = samp
            og_file_samp_assoc.save()

            assoc1 = ProcessorJobOriginalFileAssociation()
            assoc1.original_file = og_file
            assoc1.processor_job = pj
            assoc1.save()

            send_job(
                ProcessorPipeline[pj.pipeline_applied], job=pj, is_dispatch=True,
            )
            try:
                wait_for_job(pj, ProcessorJob, timezone.now())
            except Exception:
                pass

            if not pj.success:
                raise AssertionError(f"Processor job failed with reason {pj.failure_reason}")

            # Make sure that the output file agrees with our reference file

            def squish_duplicates(data: pd.DataFrame) -> pd.DataFrame:
                return data.groupby(data.index, sort=False).mean()

            # This file has already been gene-converted with the correct version
            # of the transcriptome index. We do this because the transcript IDs
            # from Ensembl change over time but the gene IDs themselves do not,
            # so we use the gene IDs as a common point of comparison.
            ref_filename = "/home/user/data_store/reference/ERR1562482_gene_converted_quant.sf"
            ref = pd.read_csv(ref_filename, delimiter="\t", index_col=0)
            ref_TPM = squish_duplicates(pd.DataFrame({"reference": ref["TPM"]}))
            ref_NumReads = squish_duplicates(pd.DataFrame({"reference": ref["NumReads"]}))

            transcript_to_gene_ids = pd.read_csv(
                "/home/user/data_store/TRANSCRIPTOME_INDEX/CAENORHABDITIS_ELEGANS/short/genes_to_transcripts.txt",
                sep="\t",
                index_col=1,
                names=["Gene"],
            )

            def convert_genes(data: pd.DataFrame) -> pd.DataFrame:
                # Map transcript IDs to gene IDs. We do this because transcript
                # IDs are not stable from version to version in Ensembl but gene
                # IDs are.
                data.index = data.index.to_series().map(
                    lambda x: transcript_to_gene_ids.loc[x, "Gene"]
                )
                return data

            output_filename = ComputedFile.objects.filter(filename="quant.sf")[0].absolute_file_path
            out = convert_genes(pd.read_csv(output_filename, delimiter="\t", index_col=0))
            out_TPM = squish_duplicates(pd.DataFrame({"actual": out["TPM"]}))
            out_NumReads = squish_duplicates(pd.DataFrame({"actual": out["NumReads"]}))

            # Make sure that there is a lot of gene overlap between the
            # reference and the output files. If this fails, it might not be the
            # end of the world but it is a good thing to know about.
            self.assertGreater(
                len(set(ref_TPM.index) & set(out_TPM.index)),
                0.95 * min(len(ref_TPM), len(out_TPM)),
            )

            (tpm_rho, _) = scipy.stats.spearmanr(ref_TPM.join(out_TPM, how="inner"))
            self.assertGreater(tpm_rho, 0.99)
            (num_reads_rho, _) = scipy.stats.spearmanr(ref_NumReads.join(out_NumReads, how="inner"))
            self.assertGreater(num_reads_rho, 0.99)


class SraRedownloadingTestCase(EndToEndTestCase):
    @tag("slow")
    @tag("salmon")
    @skip("This test is timing out I think.")
    def test_sra_redownloading(self):
        """Survey, download, then process an experiment we know is SRA."""
        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            survey_job = surveyor.survey_experiment("SRP040623", "SRA")

            self.assertTrue(survey_job.success)

            # This experiment has 4 samples that each need a downloader job.
            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 4)

            # We want one ProcessorJob to fail because it doesn't have
            # the file it was expecting, so we need to wait until one
            # DownloaderJob finishes, delete a file that is
            # downloaded, and then not delete any more.
            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            file_deleted = False
            for downloader_job in downloader_jobs:
                # We want to try and delete the file as quickly as
                # possible, so pass a short loop time and let the waiting
                # loop spin really fast so we lose as little time as
                # possible.
                downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time, 0.1)
                self.assertTrue(downloader_job.success)
                if not file_deleted:
                    for original_file in OriginalFile.objects.filter(is_downloaded=True):
                        if not original_file.is_archive:
                            original_file.delete_local_file()
                            file_deleted = True

                            # And then to make sure that we can handle
                            # cases where the downloader job is missing:
                            downloader_job.delete()
                            break

            # There's a chance that the processor job with a missing
            # file is aborted before the last downloader job
            # completes, therefore just check that there's at least 3
            # processor jobs.
            processor_jobs = ProcessorJob.objects.all()
            self.assertGreater(processor_jobs.count(), 2)

            doomed_processor_job = original_file.processor_jobs.all()[0]
            logger.info(
                "Waiting on processor Nomad job %s to fail because it realized it is missing a file.",
                doomed_processor_job.nomad_job_id,
            )

            start_time = timezone.now()
            wait_for_job(doomed_processor_job, ProcessorJob, start_time)

            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should
            # now be 5, but we also deleted on on purpose so there's 4.
            downloader_jobs = DownloaderJob.objects.all().order_by("-id")
            self.assertEqual(downloader_jobs.count(), 4)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info("Waiting on downloader Nomad job %s", recreated_job.nomad_job_id)
            recreated_job = wait_for_job(recreated_job, DownloaderJob, start_time)
            self.assertTrue(recreated_job.success)

            # Once the Downloader job succeeds, it should create one
            # and only one processor job, which the total goes back up to 4:
            self.assertEqual(ProcessorJob.objects.all().count(), 4)

            # And finally we can make sure that all of the processor
            # jobs got started correctly, including the one that got
            # recreated. However in order to save time when running
            # tests, we don't actually want to run the full salmon
            # processor. Therefore we don't have the transcriptome
            # index that is needed for this organism so the jobs will
            # fail, but that failure happens past the point that we're
            # testing.
            # So we're gonna check for the correct failure_reason.
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            good_failure_reason = "Missing transcriptome index."
            successful_processor_jobs = []
            for processor_job in processor_jobs:
                # One of the two calls to wait_for_job will fail
                # because the job is going to abort when it
                # finds that the file it wants to process is missing.
                try:
                    processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                    if not processor_job.success and processor_job.failure_reason.startswith(
                        good_failure_reason
                    ):
                        successful_processor_jobs.append(processor_job)
                except Exception:
                    pass

            self.assertEqual(len(successful_processor_jobs), 4)


class EnaFallbackTestCase(EndToEndTestCase):
    @tag("slow")
    @tag("salmon")
    @vcr.use_cassette(
        os.path.join(CASSETTES_DIR, "surveyor.test_end_to_end.unmated_reads.yaml"),
        ignore_hosts=["nomad"],
    )
    def test_unmated_reads(self):
        """Survey, download, then process a sample we know is SRA and has unmated reads.

        This test uses VCR to remove the dependence upon NCBI's
        servers, but the downloader job hits ENA's FTP and aspera
        servers. Unfortunately there's not much that can be done to
        avoid that behavior from here because the downloader jobs
        always check ENA's FTP server to see if the file has an
        unmated read. For now we'll just have to be content with the
        fact that NCBI going down won't affect this test.
        """
        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            # Survey just a single run to make things faster!
            # This sample has unmated reads!
            survey_job = surveyor.survey_experiment("SRR1603661", "SRA")

            self.assertTrue(survey_job.success)

            # Let's give the downloader a little bit to get started
            # and to update the OriginalFiles' source_urls.
            time.sleep(60)

            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 1)
            downloader_job = downloader_jobs.first()

            self.assertIsNotNone(downloader_job.start_time)

            for original_file in downloader_job.original_files.all():
                self.assertIn(".fastq.gz", original_file.source_url)

            # The downloader job will take a while to complete. Let's not wait.
            print(downloader_job.kill_nomad_job())


# This test uses the special tag "manual" because it should only be run from the "test_survey.sh" script
class SurveyTestCase(EndToEndTestCase):
    @tag("manual")
    def test_survey(self):
        """Survey the given sample"""

        # Clear out pre-existing work dirs so there's no conflicts:
        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "False")
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            survey_job = surveyor.survey_experiment(
                get_env_variable("ACCESSION"), get_env_variable("SURVEYOR")
            )

            self.assertTrue(survey_job.success)

            downloader_jobs = DownloaderJob.objects.all()
            self.assertGreater(downloader_jobs.count(), 0)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            for downloader_job in downloader_jobs:
                downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time)
                self.assertTrue(downloader_job.success)

            processor_jobs = ProcessorJob.objects.all().exclude(
                abort=True
            )  # exclude aborted processor jobs
            self.assertGreater(processor_jobs.count(), 0)

            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            start_time = timezone.now()
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                if not processor_job.success:
                    logger.error(processor_job.failure_reason)
                self.assertTrue(processor_job.success)
