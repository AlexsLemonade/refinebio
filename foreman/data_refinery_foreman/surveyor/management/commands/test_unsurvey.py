import json
import datetime
from unittest.mock import Mock, patch, call
from django.test import TransactionTestCase
from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism,
    Experiment,
    ExperimentSampleAssociation,
    Sample,
    OriginalFile,
    OriginalFileSampleAssociation
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment
from data_refinery_foreman.surveyor.geo import GeoSurveyor

class SurveyTestCase(TransactionTestCase):
    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_geo_survey_microarray(self, mock_send_task):
        """Test that the unsurveyor works correctly.

        This includes not deleting samples which also belong to other
        experiments. Therefore we survey a superseries and one of its
        sub-experiments, then delete the superseries to make sure the
        sub-experiment wasn't touched.

        We mock out the send_job function so that we don't actually
        process these. The unsurveyor code related to ComputedFile,
        ComputationalResult, and ProcessorJobs won't be tested by
        this, but it's been functionally tested.
        """
        superseries_accession = "GSE59795"
        sub_experiment_accession = "GSE46580"

        # Survey the superseries.
        survey_job = SurveyJob(source_type="GEO")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value=superseries_accession)
        key_value_pair.save()

        geo_surveyor = GeoSurveyor(survey_job)
        geo_surveyor.survey()

        # Survey the sub-experiment
        survey_job = SurveyJob(source_type="GEO")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value=sub_experiment_accession)
        key_value_pair.save()

        geo_surveyor = GeoSurveyor(survey_job)
        geo_surveyor.survey()

        # Purge the superseries
        purge_experiment(superseries_accession)

        # Make sure the subexperiment samples weren't affected.
        experiment = Experiment.objects.filter(accession_code=sub_experiment_accession)[0]
        experiment_sample_assocs = ExperimentSampleAssociation.objects.filter(experiment=experiment)
        samples = Sample.objects.filter(id__in=experiment_sample_assocs.values('sample_id'))
        self.assertEqual(samples.count(), 4)

        # Make sure sub-experiment original files weren't affected.
        og_file_sample_assocs = OriginalFileSampleAssociation.objects.filter(sample_id__in=samples.values('id'))
        original_files = OriginalFile.objects.filter(id__in=og_file_sample_assocs.values('original_file_id'))
        self.assertEqual(original_files.count(), 4)

        # Make sure the superseries samples are gone.
        self.assertEqual(Sample.objects.count(), 4)
        self.assertEqual(OriginalFile.objects.count(), 4)
