import zipfile
from datetime import timedelta
from os import path
from time import sleep
from unittest import TestCase

from django.test import tag
from django.utils import timezone

import pandas as pd
import pyrefinebio
import scipy.stats

from data_refinery_common.enums import ProcessorEnum
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    DownloaderJob,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_foreman.foreman.management.commands.create_compendia import create_compendia
from data_refinery_foreman.foreman.management.commands.create_quantpendia import create_quantpendia
from data_refinery_foreman.foreman.management.commands.run_tximport import (
    run_tximport_for_all_eligible_experiments,
)
from data_refinery_foreman.surveyor.management.commands.dispatch_qn_jobs import (
    dispatch_qn_job_if_eligible,
)
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment

SMASHER_SAMPLES = ["GSM1487313", "SRR332914"]
SMASHER_EXPERIMENTS = ["GSE1487313", "SRP332914"]

MICROARRAY_ACCESSION_CODES = [
    # TODO(ark): enable after ArrayExpress migration is done.
    # "E-TABM-496",  # 39 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE94793",  # 24 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE80822",  # 12 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE96849",  # 68 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE41094",  # 18 samples of SACCHAROMYCES_CEREVISIAE submitter processed data
]

RNA_SEQ_ACCESSION_CODES = [
    "SRP047410",  # 26 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data, one will fail.
    "SRP094706",  # 4 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
]

# TODO: remove after array express samples are enabled.
ARRAY_EXPRESS_SAMPLES_COUNT = 39

EXPECTED_TOTAL_SAMPLES_COUNT = 191 - ARRAY_EXPRESS_SAMPLES_COUNT

EXPERIMENT_ACCESSION_CODES = MICROARRAY_ACCESSION_CODES + RNA_SEQ_ACCESSION_CODES

TEST_DATA_BUCKET = "data-refinery-test-assets"


def wait_for_job(job) -> bool:
    """Waits for a job and all of its retries."""
    job.refresh_from_db()
    is_done = False

    while not is_done:

        if job.end_time is None and job.success is None:
            print(f"Polling {type(job).__name__}s. Currently waiting for job id: {job.id}")
            sleep(20)
            job.refresh_from_db()
        elif job.retried and job.retried_job:
            job = job.retried_job
        elif job.success:
            return True
        else:
            print(f"{type(job).__name__} {job.id} failed!")
            return False

    return False


def prepare_computed_files():
    # MICROARRAY TECH
    experiment = Experiment()
    experiment.accession_code = "GSE1487313"
    experiment.num_processed_samples = 1
    experiment.save()

    result = ComputationalResult()
    result.save()

    gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS", taxonomy_id=1001)

    sample = Sample()
    sample.accession_code = "GSM1487313"
    sample.title = "GSM1487313"
    sample.organism = gallus_gallus
    sample.technology = "MICROARRAY"
    sample.is_processed = True
    sample.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

    esa = ExperimentSampleAssociation()
    esa.experiment = experiment
    esa.sample = sample
    esa.save()

    computed_file = ComputedFile()
    computed_file.filename = "GSM1487313_liver.PCL"
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.is_smashable = True
    computed_file.s3_key = "GSM1487313_liver.PCL"
    computed_file.s3_bucket = TEST_DATA_BUCKET
    computed_file.save()

    assoc = SampleComputedFileAssociation()
    assoc.sample = sample
    assoc.computed_file = computed_file
    assoc.save()

    # RNASEQ TECH
    experiment2 = Experiment()
    experiment2.accession_code = "SRP332914"
    experiment2.num_processed_samples = 1
    experiment2.save()

    result2 = ComputationalResult()
    result2.save()

    sample2 = Sample()
    sample2.accession_code = "SRR332914"
    sample2.title = "SRR332914"
    sample2.organism = gallus_gallus
    sample2.technology = "RNA-SEQ"
    sample2.is_processed = True
    sample2.save()

    sra2 = SampleResultAssociation()
    sra2.sample = sample2
    sra2.result = result2
    sra2.save()

    esa2 = ExperimentSampleAssociation()
    esa2.experiment = experiment2
    esa2.sample = sample2
    esa2.save()

    computed_file2 = ComputedFile()
    computed_file2.filename = "SRP149598_gene_lengthScaledTPM.tsv"
    computed_file2.result = result2
    computed_file2.size_in_bytes = 234
    computed_file2.is_smashable = True
    computed_file2.s3_key = "SRP149598_gene_lengthScaledTPM.tsv"
    computed_file2.s3_bucket = TEST_DATA_BUCKET
    computed_file2.save()

    assoc2 = SampleComputedFileAssociation()
    assoc2.sample = sample2
    assoc2.computed_file = computed_file2
    assoc2.save()


# Use unittest TestCase instead of django TestCase to avoid the test
# being done in a transaction.
class SmasherEndToEndTestCase(TestCase):
    """Test only the smasher using precomuted samples."""

    @tag("end_to_end")
    def test_smasher_job(self):
        for accession_code in SMASHER_EXPERIMENTS:
            purge_experiment(accession_code)

        prepare_computed_files()

        # The API sometimes takes a bit to come back up.
        start_time = timezone.now()
        while True:
            try:
                pyrefinebio.create_token(agree_to_terms=True, save_token=False)
                break
            except pyrefinebio.ServerError:
                if timezone.now() - start_time > timedelta(minutes=15):
                    raise AssertionError("Server not up after 15 minutes")
                else:
                    sleep(30)

        dataset_path = "end_to_end_test_dataset"
        pyrefinebio.download_dataset(
            dataset_path,
            "testendtoend@example.com",
            dataset_dict={"GSE1487313": ["GSM1487313"], "SRP332914": ["SRR332914"]},
            timeout=timedelta(minutes=15),
            prompt=False,
        )
        self.assertTrue(path.exists(dataset_path))


# Use unittest TestCase instead of django TestCase to avoid the test
# being done in a transaction.
class FullFlowEndToEndTestCase(TestCase):
    """In order to parallelize the jobs as much as possible, everything is
    done in one big ol' function.

    This includes, in order:
      * Purging previously surveyed experiments.
      * Surveying experiments which will trigger:
        * An array_express downloader job triggering a affymetrix processor job.
        * 3 GEO downloader job triggering affymetrix processor jobs, for a total
          greater than 100 samples on a single platform.
        * A GEO downloader job triggering a NO_OP processor job.
      * Creating a transcriptome index for SACCHAROMYCES_CEREVISIAE.
      * Surveying experiments which will trigger:
        * A SRA downloader job triggering a salmon processor job for 26 samples.
          (One of which will fail, requiring a run_tximport job.)
        * A SRA downloader job triggering a salmon processor job for 4 samples.
          (This should be fully processed without intervention.)
      * Running tximport to process the experiment that had one bad sample.
      * Creating a QN Target for SACCHAROMYCES_CEREVISIAE.
      * Creating a Compendium for SACCHAROMYCES_CEREVISIAE.
      * Creating a Quantpendium for SACCHAROMYCES_CEREVISIAE.
    """

    @tag("end_to_end")
    def test_all_the_things(self):
        for accession_code in EXPERIMENT_ACCESSION_CODES:
            purge_experiment(accession_code)

        self.process_experiments()

        self.check_transcriptome_index()

        self.create_qn_reference()

        self.create_compendia()

    def process_experiments(self):
        survey_jobs = []
        # Kick off microarray jobs first because downloading the
        # affymetrix image takes a long time.
        for accession_code in MICROARRAY_ACCESSION_CODES:
            survey_jobs.append(queue_surveyor_for_accession(accession_code))

        # However, before we kick of RNA-Seq jobs, we have to build a transcriptome index for them.
        transcriptome_survey_job = queue_surveyor_for_accession("SACCHAROMYCES_CEREVISIAE, Ensembl")

        print(
            "First, creating transcriptome indices (and starting on Affy jobs while we're at it):"
        )
        self.assertTrue(wait_for_job(transcriptome_survey_job))

        transcriptome_downloader_jobs = DownloaderJob.objects.filter(
            downloader_task="TRANSCRIPTOME_INDEX",
            created_at__gt=transcriptome_survey_job.created_at,
        )

        for downloader_job in transcriptome_downloader_jobs:
            self.assertTrue(wait_for_job(downloader_job))

        transcriptome_processor_jobs = ProcessorJob.objects.filter(
            pipeline_applied__startswith="TRANSCRIPTOME_INDEX",
            created_at__gt=transcriptome_survey_job.created_at,
        )

        for processor_job in transcriptome_processor_jobs:
            self.assertTrue(wait_for_job(processor_job))

        for accession_code in RNA_SEQ_ACCESSION_CODES:
            survey_jobs.append(queue_surveyor_for_accession(accession_code))

        print("Next, processing all the raw data.")
        for survey_job in survey_jobs:
            self.assertTrue(wait_for_job(survey_job))

        # We need to ignore a lot of samples in staging that aren't
        # related to these experiments.
        # sample_id_list = Experiment.objects.filter(accession_code__in=EXPERIMENT_ACCESSION_CODES)
        sample_id_list = Sample.objects.filter(
            experiment__accession_code__in=EXPERIMENT_ACCESSION_CODES
        ).values_list("id")

        self.assertEqual(
            Sample.objects.filter(id__in=sample_id_list).count(), EXPECTED_TOTAL_SAMPLES_COUNT
        )

        samples = []
        for accession_code in EXPERIMENT_ACCESSION_CODES:
            experiment = Experiment.objects.get(accession_code=accession_code)
            samples.extend(list(experiment.samples.all()))

        downloader_jobs = []
        for sample in samples:
            greatest_job_id = -1
            last_job = None

            # There should be only one, but if it got retried for some
            # reason we want the latest one.
            for job in sample.get_downloader_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            downloader_jobs.append(last_job)

        for downloader_job in downloader_jobs:
            self.assertTrue(wait_for_job(downloader_job))

        processor_jobs = []
        for sample in samples:
            # This sample fails for good reason, so don't expect it to pass.
            if sample.accession_code == "SRR1583739":
                continue

            greatest_job_id = -1
            last_job = None

            for job in sample.get_processor_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            processor_jobs.append(last_job)

        for processor_job in processor_jobs:
            self.assertTrue(wait_for_job(processor_job))

        # Because SRR1583739 fails, the 26 samples from SRP047410 won't be processed
        EXPECTED_SAMPLE_FAILURES = 26
        self.assertEqual(
            Sample.processed_objects.filter(id__in=sample_id_list).count(),
            EXPECTED_TOTAL_SAMPLES_COUNT - EXPECTED_SAMPLE_FAILURES,
        )

        print("Finally, need to run tximport to finish an experiment with one bad sample.")
        tximport_jobs = run_tximport_for_all_eligible_experiments(dispatch_jobs=False)
        self.assertEqual(len(tximport_jobs), 1)

        self.assertTrue(wait_for_job(tximport_jobs[0]))

        # This is the full total of jobs minus one.
        self.assertEqual(
            Sample.processed_objects.filter(id__in=sample_id_list).count(),
            EXPECTED_TOTAL_SAMPLES_COUNT - 1,
        )

    def check_transcriptome_index(self):
        # Make sure that a processed file using our new transcriptome index is
        # similar enough to a reference file
        print(
            "Now we are going to verify that the outputs of salmon look okay with"
            "the transcriptome index we generated."
        )

        sample = Sample.objects.get(accession_code="SRR5085168")

        # First, sanity check that the genome build hasn't changed. If it
        # changes, all bets are off when comparing quant.sf files.
        self.assertEqual(
            sample.results.all()
            .filter(processor__name=ProcessorEnum.SALMON_QUANT.value["name"])
            .first()
            .organism_index.assembly_name,
            "R64-1-1",
        )

        # We are now going to download the `quant.sf` file by making a
        # quant-only smasher job. We need to do this because the results bucket
        # was locked down, so we can't access it without agreeing to the terms
        # first.
        pyrefinebio.create_token(agree_to_terms=True, save_token=False)

        dataset_path = "end_to_end_quant_dataset.zip"
        pyrefinebio.download_dataset(
            dataset_path,
            "testendtoend@example.com",
            dataset_dict={"SRP094706": ["SRR5085168"]},
            quant_sf_only=True,
            aggregation="EXPERIMENT",
            timeout=timedelta(minutes=15),
            prompt=False,
        )
        self.assertTrue(path.exists(dataset_path))

        def squish_duplicates(data: pd.DataFrame) -> pd.DataFrame:
            return data.groupby(data.index, sort=False).mean()

        ref_filename = "/home/user/data_store/reference/SRR5085168_quant.sf"
        ref = pd.read_csv(ref_filename, delimiter="\t", index_col=0)
        ref_TPM = squish_duplicates(pd.DataFrame({"reference": ref["TPM"]}))
        ref_NumReads = squish_duplicates(pd.DataFrame({"reference": ref["NumReads"]}))

        with zipfile.ZipFile(dataset_path) as zf:
            with zf.open("SRP094706/SRR5085168_quant.sf", "r") as out_file:
                out = pd.read_csv(out_file, delimiter="\t", index_col=0)
        out_TPM = squish_duplicates(pd.DataFrame({"actual": out["TPM"]}))
        out_NumReads = squish_duplicates(pd.DataFrame({"actual": out["NumReads"]}))

        # Make sure that there is a lot of gene overlap between the
        # reference and the output files using a Jaccard index test.
        intersection = len(set(ref_TPM.index) & set(out_TPM.index))
        union = len(set(ref_TPM.index) | set(out_TPM.index))
        jaccard_index = intersection / union
        self.assertGreater(jaccard_index, 0.95)

        (tpm_rho, _) = scipy.stats.spearmanr(ref_TPM.join(out_TPM, how="inner"))
        self.assertGreater(tpm_rho, 0.99)
        (num_reads_rho, _) = scipy.stats.spearmanr(ref_NumReads.join(out_NumReads, how="inner"))
        self.assertGreater(num_reads_rho, 0.99)

    def create_qn_reference(self):
        organism = Organism.objects.get(name="SACCHAROMYCES_CEREVISIAE")

        qn_job = dispatch_qn_job_if_eligible(organism)
        self.assertIsNotNone(qn_job)

        self.assertTrue(wait_for_job(qn_job))

    def create_compendia(self):
        """Create both a compendium and a quantpendium"""
        compendia_jobs = create_compendia(None, "SACCHAROMYCES_CEREVISIAE")
        quantpendia_jobs = create_quantpendia("SACCHAROMYCES_CEREVISIAE", None)

        for job in compendia_jobs + quantpendia_jobs:
            self.assertTrue(wait_for_job(job))
