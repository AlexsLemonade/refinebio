import json
import datetime
from unittest.mock import Mock, patch, call
from django.test import TransactionTestCase
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism,
    Sample
)
from data_refinery_foreman.surveyor.geo import (
    GeoSurveyor
)


class SurveyTestCase(TransactionTestCase):

    def prep_test(self, experiment_accession):
        survey_job = SurveyJob(source_type="GEO")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value=experiment_accession)
        key_value_pair.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_geo_survey_microarray(self, mock_send_task):
        """ Run the GEO surveyor and make sure we get some files to DL!

        For an Illumina Microarray platform.
        """
        self.prep_test("GSE11915")

        geo_surveyor = GeoSurveyor(self.survey_job)
        geo_surveyor.survey()

        self.assertEqual(34, Sample.objects.all().count())

        sample_object = Sample.objects.first()
        self.assertEqual(sample_object.platform_name, "[HG-U133A] Affymetrix Human Genome U133A Array")
        self.assertEqual(sample_object.platform_accession_code, "hgu133a")
        self.assertEqual(sample_object.technology, "MICROARRAY")

        downloader_jobs = DownloaderJob.objects.all()
        self.assertEqual(46, downloader_jobs.count())

    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_geo_survey_agilent(self, mock_send_task):
        """ Run the GEO surveyor and make sure we get some files to DL!

        For an Agilent Microarray platform.
        """
        self.prep_test("GSE35186")

        geo_surveyor = GeoSurveyor(self.survey_job)
        geo_surveyor.survey()

        self.assertEqual(124, Sample.objects.all().count())

        sample_object = Sample.objects.first()
        self.assertEqual(sample_object.platform_name,
                         "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)")
        self.assertEqual(sample_object.platform_accession_code, "GPL6480")
        # We currently do not support Agilent platforms, so we can't
        # match its accession to one we know about.
        self.assertEqual(sample_object.technology, "UNKNOWN")

        downloader_jobs = DownloaderJob.objects.all()
        # There would be 124 samples + 2 metadata files. However at
        # the moment Agilent is unsupported so we don't want to queue
        # downloader jobs.
        self.assertEqual(0, downloader_jobs.count())

    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_geo_survey_rnaseq(self, mock_send_task):
        """ Run the GEO surveyor and make sure we get some files to DL!

        For an Illumina RNASeq platform.
        """
        self.prep_test("GSE99264")

        geo_surveyor = GeoSurveyor(self.survey_job)
        geo_surveyor.survey()

        self.assertEqual(7, Sample.objects.all().count())

        sample_object = Sample.objects.first()
        self.assertEqual(sample_object.platform_name, "Illumina Genome Analyzer II")
        self.assertEqual(sample_object.platform_accession_code, "Illumina Genome Analyzer II")
        self.assertEqual(sample_object.technology, "RNA-SEQ")

        downloader_jobs = DownloaderJob.objects.all()
        self.assertEqual(2, downloader_jobs.count())

    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_geo_survey_superseries(self, mock_send_task):
        """ Run the GEO surveyor and make sure we get some files to DL!

        For a Super Series.
        """
        self.prep_test("GSE54334")

        geo_surveyor = GeoSurveyor(self.survey_job)
        geo_surveyor.survey()

        self.assertEqual(7, Sample.objects.all().count())

        sample_object = Sample.objects.first()
        self.assertEqual(sample_object.platform_name, "Illumina HiSeq 2500")
        self.assertEqual(sample_object.platform_accession_code, "Illumina HiSeq 2500")
        self.assertEqual(sample_object.technology, "RNA-SEQ")

        downloader_jobs = DownloaderJob.objects.all()
        self.assertEqual(2, downloader_jobs.count())
