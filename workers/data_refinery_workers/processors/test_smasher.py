import os
from django.test import TestCase, tag
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    ComputationalResult,
    ComputedFile,
    Experiment,
    Organism,
    Sample,
    SampleResultAssociation,
    Dataset,
    ProcessorJobDatasetAssociation
)
from data_refinery_workers.processors import smasher
from data_refinery_workers.processors import utils

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SMASHER"
    pj.save()

    experiment = Experiment()
    experiment.accession_code = "GSE51081"
    experiment.save()

    result = ComputationalResult()
    result.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

    sample = Sample()
    sample.accession_code = 'GSM1237810'
    sample.organism = homo_sapiens
    sample.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

    computed_file = ComputedFile()
    computed_file.filename = "GSM1237810_T09-1084.PCL"
    computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.save()

    sample = Sample()
    sample.accession_code = 'GSM1237812'
    sample.organism = homo_sapiens
    sample.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

    computed_file = ComputedFile()
    computed_file.filename = "GSM1237812_S97-PURE.PCL"
    computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.save()

    ds = Dataset()
    ds.data = {'GSE51081': ['GSM1237810', 'GSM1237812']}
    ds.aggregate_by = 'EXPERIMENT' # [ALL or SPECIES or EXPERIMENT]
    ds.email_address = "null@null.null"
    ds.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = pj
    pjda.dataset = ds
    pjda.save()

    return pj

class SmasherTestCase(TestCase):

    @tag("smasher")
    def test_smasher(self):
        """ Main tester. """
        job = prepare_job()

        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)
        dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
        job_context_check = {}
        job_context_check["dataset"] = dataset
        job_context_check["samples"] = dataset.get_samples()
        job_context_check["experiments"] = dataset.get_experiments()
        self.assertEqual(len(job_context_check['samples']), 2)
        self.assertEqual(len(job_context_check['experiments']), 1)

        # Smoke test while we're here..
        dataset.get_samples_by_experiment()
        dataset.get_samples_by_species()
        dataset.get_aggregated_samples()

        for ag_type in ['ALL', 'EXPERIMENT', 'SPECIES']:
            dataset.aggregate_by = ag_type
            dataset.save()

            final_context = smasher.smash(job.pk)
            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context['output_file']), 0)
            self.assertEqual(final_context['dataset'].is_processed, True)

            dataset.is_processed = False
            dataset.save()

            # Cleanup
            os.remove(final_context['output_file']) 

    @tag("smasher")
    def test_get_results(self):
        """ Test our ability to collect the appropriate samples. """

        sample = Sample()
        sample.accession_code = 'GSM45588'
        sample.save()

        result = ComputationalResult()
        result.save()

        computed_file = ComputedFile()
        computed_file.filename = "oh_boy.txt"
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.save()

        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        computed_files = sample.get_result_files()
        self.assertEqual(computed_files.count(), 1)

    @tag("smasher")
    def test_sanity_imports(self):
        """ Sci imports can be tricky, make sure this works. """

        import numpy
        import scipy
        import matplotlib
        import pandas
        import sklearn
        import sympy
