from django.test import TestCase

from data_refinery_common.models import Experiment, ExperimentSampleAssociation, Sample


class ExperimentModelTestCase(TestCase):
    def tearDown(self):
        Experiment.objects.all().delete()
        Sample.objects.all().delete()

    # Test for `get_sample_metadata_fields`
    def test_get_sample_metadata_fields(self):
        experiment = Experiment()
        experiment.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.specimen_part = "Lung"
        sample.sex = "Male"
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        self.assertEqual(
            set(experiment.get_sample_metadata_fields()), set(["specimen_part", "sex"])
        )

    # Test for when no metadata fields are present
    def test_get_sample_metadata_fields_none(self):
        experiment = Experiment()
        experiment.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        self.assertEqual(experiment.get_sample_metadata_fields(), [])

    # Since 'age' is stored as a number, it can lead to errors where other fields wouldn't
    def test_get_sample_metadata_fields_age(self):
        experiment = Experiment()
        experiment.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.age = 23
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        self.assertEqual(set(experiment.get_sample_metadata_fields()), set(["age"]))
