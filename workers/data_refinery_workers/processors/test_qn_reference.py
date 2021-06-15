import os
import os.path
from io import StringIO
from typing import List

from django.core.management import call_command
from django.test import TestCase, tag

import numpy as np

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SampleComputedFileAssociation,
)
from data_refinery_common.models.organism import Organism
from data_refinery_workers.processors import qn_reference, smasher


def prepare_experiment(ids: List[int]) -> Experiment:
    (homo_sapiens, _) = Organism.objects.get_or_create(name="HOMO_SAPIENS", taxonomy_id=9606)

    experiment = Experiment()
    experiment.accession_code = "12345"
    experiment.save()
    codes = [str(i) for i in ids]

    for code in codes:
        sample = Sample()
        sample.accession_code = code
        sample.title = code
        sample.platform_accession_code = "A-MEXP-1171"
        sample.manufacturer = "SLIPPERY DICK'S DISCOUNT MICROARRAYS"
        sample.organism = homo_sapiens
        sample.technology = "MICROARRAY"
        sample.is_processed = True
        sample.save()

        cr = ComputationalResult()
        cr.save()

        computed_file = ComputedFile()
        computed_file.filename = code + ".tsv"
        computed_file.absolute_file_path = "/home/user/data_store/QN/" + code + ".tsv"
        computed_file.size_in_bytes = int(code)
        computed_file.result = cr
        computed_file.is_smashable = True
        computed_file.save()

        scfa = SampleComputedFileAssociation()
        scfa.sample = sample
        scfa.computed_file = computed_file
        scfa.save()

        exsa = ExperimentSampleAssociation()
        exsa.experiment = experiment
        exsa.sample = sample
        exsa.save()


# Calculated using excel
EXPECTED_REFERENCE_LIST = [
    -8.1385265774612,
    -6.63998543415348,
    -4.75403705778203,
    -3.76609255353842,
    -2.82302739486148,
    -1.97927983461526,
    -1.09596077865308,
    -0.437948852881293,
    -0.064011709095852,
    0.445693201406115,
    1.17059982897098,
    1.51254730833957,
    2.29756769541054,
    2.53656272830995,
    2.87634550809431,
    3.39628303134281,
    3.6771607226078,
    3.93075578415268,
    4.44528396687628,
    5.03900163488493,
    5.23997418236877,
    5.68127443617785,
    6.09602667800399,
    6.84706563853497,
    7.73280620381975,
    7.98338751295354,
    8.204330858066,
    9.35870852286352,
    10.3978144037197,
    13.9271586069496,
]


class QNRefTestCase(TestCase):
    @tag("qn")
    def test_qn_reference(self):
        # We don't have a 0.tsv
        experiment = prepare_experiment(range(1, 201))

        job = ProcessorJob()
        job.pipeline_applied = "QN_REFERENCE"
        job.save()

        dataset = Dataset()
        dataset.data = {"12345": ["1", "2", "3", "4", "5", "6"]}
        dataset.aggregate_by = "ALL"
        dataset.scale_by = "NONE"
        dataset.quantile_normalize = False  # We don't QN because we're creating the target now
        dataset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        final_context = qn_reference.create_qn_reference(job.pk)
        self.assertTrue(final_context["success"])
        self.assertTrue(os.path.exists(final_context["target_file"]))
        self.assertEqual(os.path.getsize(final_context["target_file"]), 562)

        homo_sapiens = Organism.objects.get(taxonomy_id=9606)
        target = homo_sapiens.qn_target.computedfile_set.latest()
        self.assertEqual(target.sha1, "de69d348f8b239479e2330d596c4013a7b0b2b6a")

        # Create and run a smasher job that will use the QN target we just made.
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        ds = Dataset()
        ds.data = {"12345": ["1", "2", "3", "4", "5"]}
        ds.aggregate_by = "SPECIES"
        ds.scale_by = "STANDARD"
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = True
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(final_context["success"])

        np.testing.assert_almost_equal(final_context["merged_qn"]["1"][0], -0.4379488527774811)
        np.testing.assert_almost_equal(final_context["original_merged"]["1"][0], -0.5762109)

        # Make sure that the results were created. We create 200 computed files
        # and computational results (1 for each sample) plus the one generated
        # by the QN reference processor.
        self.assertEqual(ComputedFile.objects.all().count(), 200 + 1)
        self.assertEqual(ComputationalResult.objects.all().count(), 200 + 1)
        self.assertEqual(ComputationalResultAnnotation.objects.all().count(), 1)

    @tag("qn")
    def test_build_qn_target_correct_files(self):
        # We don't have a 0.tsv
        experiment = prepare_experiment(range(1, 201))

        job_context = {
            "target_file": "/home/user/data_store/QN/output/build_qn_target_correct_files.tsv",
            "input_files": {"ALL": []},
        }

        # 6.tsv and 7.tsv have invalid data, so we'll skip them for this test
        for i in range(1, 6):
            sample = Sample.objects.get(accession_code=i)
            job_context["input_files"]["ALL"].extend(
                (f, sample) for f in sample.computed_files.all()
            )

        target_dir = os.path.dirname(job_context["target_file"])
        os.makedirs(target_dir, exist_ok=True)

        job_context = qn_reference._build_qn_target(job_context)

        self.assertEquals(job_context["geneset"], list(set(map(str, range(1, 31)))))
        self.assertEquals(job_context["num_valid_inputs"], 5)

        sum_frame = job_context["sum_frame"]
        sum_frame.index = sum_frame.index.to_series().map(int)
        sum_list = list(sum_frame["sum"])

        for i in range(0, 30):
            self.assertAlmostEqual(sum_list[i], EXPECTED_REFERENCE_LIST[i], delta=0.0001)

    @tag("qn")
    def test_build_qn_invalid_files(self):
        # We don't have a 0.tsv
        experiment = prepare_experiment(range(1, 201))

        job_context = {
            "target_file": "/home/user/data_store/QN/output/build_qn_target_invalid_files.tsv",
            "input_files": {"ALL": []},
        }

        # 6.tsv has an invalid gene (31)
        # 7.tsv has some NA values
        for i in range(1, 8):
            sample = Sample.objects.get(accession_code=i)
            job_context["input_files"]["ALL"].extend(
                (f, sample) for f in sample.computed_files.all()
            )

        target_dir = os.path.dirname(job_context["target_file"])
        os.makedirs(target_dir, exist_ok=True)

        job_context = qn_reference._build_qn_target(job_context)

        self.assertEquals(job_context["geneset"], list(set(map(str, range(1, 31)))))
        self.assertEquals(job_context["num_valid_inputs"], 5)

        sum_frame = job_context["sum_frame"]
        sum_frame.index = sum_frame.index.to_series().map(int)
        sum_list = list(sum_frame["sum"])

        for i in range(0, 30):
            self.assertAlmostEqual(sum_list[i], EXPECTED_REFERENCE_LIST[i], delta=0.0001)

    @tag("qn")
    def test_qn_management_command(self):
        """Test that the management command fires off and then does not create
        a job for an organism that does not have enough samples on the same
        platform."""

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        homo_sapiens.save()

        # We don't have a 0.tsv
        experiment = prepare_experiment(range(1, 7))

        out = StringIO()
        try:
            call_command("create_qn_target", organism="homo_sapiens", min=1, stdout=out)
        except SystemExit as e:  # this is okay!
            pass

        stdout = out.getvalue()
        self.assertFalse("Target file" in stdout)

        # There's not enough samples available in this scenario so we
        # shouldn't have even made a processor job.
        self.assertEqual(ProcessorJob.objects.count(), 0)
