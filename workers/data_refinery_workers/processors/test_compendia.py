import json
import os
import zipfile
from typing import Dict

from django.test import TransactionTestCase, tag

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    SurveyJob,
)
from data_refinery_workers.processors import create_compendia


def create_sample_for_experiment(sample_info: Dict, experiment: Experiment) -> Sample:
    result = ComputationalResult()
    result.save()

    sample = Sample()
    sample.accession_code = sample_info["accession_code"]
    sample.title = sample_info.get("title", None) or sample_info["accession_code"]
    sample.organism = sample_info["organism"]
    sample.technology = sample_info["technology"]
    sample.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

    esa = ExperimentSampleAssociation()
    esa.experiment = experiment
    esa.sample = sample
    esa.save()

    if sample_info.get("filename") is not None:
        computed_file = ComputedFile()
        computed_file.filename = sample_info["filename"]
        computed_file.absolute_file_path = sample_info["data_dir"] + sample_info["filename"]
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

    return sample


class CompendiaTestCase(TransactionTestCase):
    @tag("compendia")
    def test_create_compendia(self):
        DATA_DIR = "/home/user/data_store/PCL/"

        job = ProcessorJob()
        job.pipeline_applied = ProcessorPipeline.CREATE_COMPENDIA.value
        job.save()

        gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS", taxonomy_id=1001)

        # MICROARRAY TECH
        (experiment, _) = Experiment.objects.get_or_create(accession_code="GSE1487313")
        experiment.accession_code = "GSE1487313"
        experiment.save()

        create_sample_for_experiment(
            {
                "organism": gallus_gallus,
                "accession_code": "GSM1487313",
                "technology": "MICROARRAY",
                "filename": "GSM1487313_liver.PCL",
                "data_dir": DATA_DIR,
            },
            experiment,
        )

        # Missing sample that will be filtered
        create_sample_for_experiment(
            {
                "organism": gallus_gallus,
                "accession_code": "GSM1487222",
                "title": "this sample will be filtered",
                "technology": "MICROARRAY",
                "filename": "GSM1487222_empty.PCL",
                "data_dir": DATA_DIR,
            },
            experiment,
        )

        # RNASEQ TECH
        experiment2 = Experiment()
        experiment2.accession_code = "SRS332914"
        experiment2.save()

        create_sample_for_experiment(
            {
                "organism": gallus_gallus,
                "accession_code": "SRS332914",
                "technology": "RNA-SEQ",
                "filename": "SRP149598_gene_lengthScaledTPM.tsv",
                "data_dir": DATA_DIR,
            },
            experiment,
        )

        dset = Dataset()
        dset.data = {
            "GSE1487313": ["GSM1487313", "GSM1487222"],
            "SRX332914": ["SRS332914"],
        }
        dset.scale_by = "NONE"
        dset.aggregate_by = "SPECIES"
        dset.svd_algorithm = "ARPACK"
        dset.quantile_normalize = True
        dset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dset
        pjda.save()

        final_context = create_compendia.create_compendia(job.id)

        # Because one of the samples is filtered out, there will be too few
        # remaining samples to smash together, so we expect this job to fail.
        job.refresh_from_db()
        self.assertFalse(job.success)
        self.assertIsNotNone(job.failure_reason)
        self.assertIn("k must be between 1 and min(A.shape)", job.failure_reason)

        # check that sample with no computed file was skipped
        self.assertTrue("GSM1487222" in final_context["filtered_samples"])
        self.assertEqual(
            final_context["filtered_samples"]["GSM1487222"]["experiment_accession_code"],
            "GSE1487313",
        )

    @tag("compendia")
    def test_create_compendia_danio(self):
        job = ProcessorJob()
        job.pipeline_applied = ProcessorPipeline.CREATE_COMPENDIA.value
        job.save()

        # MICROARRAY TECH
        experiment = Experiment()
        experiment.accession_code = "GSE1234"
        experiment.save()

        result = ComputationalResult()
        result.save()

        qn_target = ComputedFile()
        qn_target.filename = "danio_target.tsv"
        qn_target.absolute_file_path = "/home/user/data_store/QN/danio_target.tsv"
        qn_target.is_qn_target = True
        qn_target.size_in_bytes = "12345"
        qn_target.sha1 = "aabbccddeeff"
        qn_target.result = result
        qn_target.save()

        danio_rerio = Organism(name="DANIO_RERIO", taxonomy_id=1, qn_target=result)
        danio_rerio.save()

        cra = ComputationalResultAnnotation()
        cra.data = {}
        cra.data["organism_id"] = danio_rerio.id
        cra.data["is_qn"] = True
        cra.result = result
        cra.save()

        result = ComputationalResult()
        result.save()

        micros = []
        for file in os.listdir("/home/user/data_store/raw/TEST/MICROARRAY/"):

            if "microarray.txt" in file:
                continue

            create_sample_for_experiment(
                {
                    "organism": danio_rerio,
                    "accession_code": file,
                    "technology": "MICROARRAY",
                    "filename": file,
                    "data_dir": "/home/user/data_store/raw/TEST/MICROARRAY/",
                },
                experiment,
            )

            micros.append(file)

        experiment = Experiment()
        experiment.accession_code = "GSE5678"
        experiment.save()

        result = ComputationalResult()
        result.save()
        rnas = []
        for file in os.listdir("/home/user/data_store/raw/TEST/RNASEQ/"):

            if "rnaseq.txt" in file:
                continue

            create_sample_for_experiment(
                {
                    "organism": danio_rerio,
                    "accession_code": file,
                    "technology": "RNA-SEQ",
                    "filename": file,
                    "data_dir": "/home/user/data_store/raw/TEST/RNASEQ/",
                },
                experiment,
            )

            rnas.append(file)

        # Missing sample that will be filtered
        sample = create_sample_for_experiment(
            {
                "organism": danio_rerio,
                "accession_code": "GSM1487222",
                "title": "this sample will be filtered",
                "technology": "RNA-SEQ",
                "filename": None,
            },
            experiment,
        )
        rnas.append(sample.accession_code)

        dset = Dataset()
        dset.data = {"GSE1234": micros, "GSE5678": rnas}
        dset.scale_by = "NONE"
        dset.aggregate_by = "SPECIES"
        dset.svd_algorithm = "ARPACK"
        dset.quantile_normalize = True
        dset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dset
        pjda.save()

        final_context = create_compendia.create_compendia(job.id)

        # Verify result
        self.assertEqual(final_context["compendium_result"].result.computedfile_set.count(), 1)
        for file in final_context["compendium_result"].result.computedfile_set.all():
            self.assertTrue(os.path.exists(file.absolute_file_path))

        # test compendium_result
        self.assertEqual(final_context["compendium_result"].svd_algorithm, "ARPACK")
        self.assertEqual(
            final_context["compendium_result"].primary_organism.name,
            final_context["organism_name"],
        )
        self.assertEqual(final_context["compendium_result"].primary_organism.name, "DANIO_RERIO")
        self.assertEqual(final_context["compendium_result"].organisms.count(), 1)

        # check that sample with no computed file was skipped
        self.assertTrue("GSM1487222" in final_context["filtered_samples"])
        self.assertEqual(
            final_context["filtered_samples"]["GSM1487222"]["experiment_accession_code"], "GSE5678",
        )

        # It's maybe not worth asserting this until we're sure the behavior is correct
        # self.assertEqual(final_context['merged_qn'].shape, (9045, 830))

        zf = zipfile.ZipFile(
            final_context["compendium_result"].result.computedfile_set.first().absolute_file_path
        )
        with zf.open("aggregated_metadata.json") as f:
            metadata = json.load(f)

            self.assertFalse(metadata.get("quant_sf_only"))
            # 420 microarray + 419 RNA seq (+1 that should be filtered)
            self.assertEqual(metadata.get("num_samples"), 839)
            self.assertEqual(metadata.get("num_experiments"), 2)

            # Make sure the data were quantile normalized
            self.assertTrue(metadata.get("quantile_normalized"))
