import datetime
import json
from unittest.mock import Mock, call, patch

from django.test import TransactionTestCase

from data_refinery_common.models import (
    ComputationalResult,
    Experiment,
    ExperimentResultAssociation,
    ExperimentSampleAssociation,
    Organism,
    Processor,
    Sample,
    SampleResultAssociation,
)
from data_refinery_foreman.foreman.management.commands.assoc_experiment_results import (
    make_experiment_result_associations,
)
from data_refinery_foreman.surveyor.geo import GeoSurveyor


class SurveyTestCase(TransactionTestCase):
    def test_make_experiment_result_associations(self):
        """Tests that the correct associations are made.

        The situation we're setting up is basically this:
          * tximport has been run for an experiment.
          * It made associations between the samples in
            the experiment and the ComputationalResult.
          * It didn't make associations between the
            experiment itself and the ComputationalResult.
          * There is a second experiment that hasn't had
            tximport run but shares a sample with the
            other experiment.
          * This second experiment has a sample which has
            not yet had tximport run on it.

        And what we're going to test for is:
          * An association is created between the tximport
            result and the first experiment.
          * An association is NOT created between the
            tximport result and the second experiment.
        """
        # Get an organism to set on samples:
        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        # Create the tximport processor and result:
        processor = Processor()
        processor.name = "Tximport"
        processor.version = "v9.9.9"
        processor.docker_image = "dr_salmon"
        processor.environment = '{"some": "environment"}'
        processor.save()

        result = ComputationalResult()
        result.commands.append("tximport invocation")
        result.is_ccdl = True
        result.processor = processor
        result.save()

        # Create the first experiment and it's samples:
        processed_experiment = Experiment()
        processed_experiment.accession_code = "SRP12345"
        processed_experiment.save()

        processed_sample_one = Sample()
        processed_sample_one.accession_code = "SRX12345"
        processed_sample_one.title = "SRX12345"
        processed_sample_one.organism = homo_sapiens
        processed_sample_one.save()

        sra = SampleResultAssociation()
        sra.sample = processed_sample_one
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = processed_experiment
        esa.sample = processed_sample_one
        esa.save()

        processed_sample_two = Sample()
        processed_sample_two.accession_code = "SRX12346"
        processed_sample_two.title = "SRX12346"
        processed_sample_two.organism = homo_sapiens
        processed_sample_two.save()

        sra = SampleResultAssociation()
        sra.sample = processed_sample_two
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = processed_experiment
        esa.sample = processed_sample_two
        esa.save()

        # Create the second experiment and it's additional sample.
        unprocessed_experiment = Experiment()
        unprocessed_experiment.accession_code = "SRP6789"
        unprocessed_experiment.save()

        unprocessed_sample = Sample()
        unprocessed_sample.accession_code = "SRX6789"
        unprocessed_sample.title = "SRX6789"
        unprocessed_sample.organism = homo_sapiens
        unprocessed_sample.save()

        sra = SampleResultAssociation()
        sra.sample = unprocessed_sample
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = unprocessed_experiment
        esa.sample = unprocessed_sample
        esa.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = unprocessed_experiment
        esa.sample = processed_sample_two
        esa.save()

        # Run the function we're testing:
        make_experiment_result_associations()

        # Test that only one association was created and that it was
        # to the processed experiment:
        eras = ExperimentResultAssociation.objects.all()

        self.assertEqual(len(eras), 1)
        self.assertEqual(eras.first().experiment, processed_experiment)
