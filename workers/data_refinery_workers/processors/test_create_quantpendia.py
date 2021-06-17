import json
import os
import zipfile

from django.test import TransactionTestCase, tag

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_workers.processors.create_quantpendia import create_quantpendia


class QuantpendiaTestCase(TransactionTestCase):
    @tag("compendia")
    def test_create_quantpendia(self):
        job = ProcessorJob()
        job.pipeline_applied = ProcessorPipeline.CREATE_QUANTPENDIA.value
        job.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51088"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=9606)

        sample = Sample()
        sample.accession_code = "GSM1237818"
        sample.title = "GSM1237818"
        sample.organism = homo_sapiens
        sample.technology = "RNA-SEQ"
        sample.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = experiment
        esa.sample = sample
        esa.save()

        # Add a second non-downloadable sample. This one should not be included
        # in the count of samples available in the metadata
        sample2 = Sample()
        sample2.accession_code = "GSM1237819"
        sample2.title = "GSM1237819"
        sample2.organism = homo_sapiens
        sample2.technology = "RNA-SEQ"
        sample2.save()

        esa2 = ExperimentSampleAssociation()
        esa2.experiment = experiment
        esa2.sample = sample2
        esa2.save()

        computed_file = ComputedFile()
        computed_file.s3_key = "smasher-test-quant.sf"
        computed_file.s3_bucket = "data-refinery-test-assets"
        computed_file.filename = "quant.sf"
        computed_file.absolute_file_path = "/home/user/data_store/QUANT/smasher-test-quant.sf"
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.size_in_bytes = 123123
        computed_file.sha1 = (
            "08c7ea90b66b52f7cd9d9a569717a1f5f3874967"  # this matches with the downloaded file
        )
        computed_file.save()

        computed_file = ComputedFile()
        computed_file.filename = "logquant.tsv"
        computed_file.is_smashable = True
        computed_file.size_in_bytes = 123123
        computed_file.result = result
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {"GSE51088": ["GSM1237818", "GSM1237819"]}
        ds.aggregate_by = "EXPERIMENT"
        ds.scale_by = "STANDARD"
        ds.email_address = "null@derp.com"
        ds.quant_sf_only = True  # Make the dataset include quant.sf files only
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = create_quantpendia(job.id)

        self.assertTrue(
            os.path.exists(final_context["output_dir"] + "/GSE51088/GSM1237818_quant.sf")
        )
        self.assertTrue(os.path.exists(final_context["output_dir"] + "/README.md"))
        self.assertTrue(os.path.exists(final_context["output_dir"] + "/LICENSE.TXT"))
        self.assertTrue(os.path.exists(final_context["output_dir"] + "/aggregated_metadata.json"))

        # test that archive exists
        quantpendia_file = ComputedFile.objects.filter(
            is_compendia=True, quant_sf_only=True
        ).latest()
        self.assertTrue(os.path.exists(quantpendia_file.absolute_file_path))

        zf = zipfile.ZipFile(quantpendia_file.absolute_file_path)
        with zf.open("aggregated_metadata.json") as f:
            metadata = json.load(f)

            self.assertTrue(metadata.get("quant_sf_only"))
            self.assertEqual(metadata.get("num_samples"), 1)
            self.assertEqual(metadata.get("num_experiments"), 1)

            # Make sure the data were not quantile normalized
            self.assertFalse(metadata.get("quantile_normalized"))
