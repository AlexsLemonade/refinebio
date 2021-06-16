import json
import os
import zipfile

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


class CompendiaTestCase(TransactionTestCase):
    @tag("compendia")
    def test_sanity(self):
        print("Hey!")

    @tag("compendia")
    def test_create_compendia(self):
        job = ProcessorJob()
        job.pipeline_applied = ProcessorPipeline.CREATE_COMPENDIA.value
        job.save()

        # MICROARRAY TECH
        experiment = Experiment()
        experiment.accession_code = "GSE1487313"
        experiment.save()

        result = ComputationalResult()
        result.save()

        gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "GSM1487313"
        sample.title = "GSM1487313"
        sample.organism = gallus_gallus
        sample.technology = "MICROARRAY"
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
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        # Missing sample that will be filtered
        sample = Sample()
        sample.accession_code = "GSM1487222"
        sample.title = "this sample will be filtered"
        sample.organism = gallus_gallus
        sample.technology = "MICROARRAY"
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
        computed_file.filename = "GSM1487222_empty.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/doesnt_exists.PCL"
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        # RNASEQ TECH
        experiment2 = Experiment()
        experiment2.accession_code = "SRS332914"
        experiment2.save()

        result2 = ComputationalResult()
        result2.save()

        sample2 = Sample()
        sample2.accession_code = "SRS332914"
        sample2.title = "SRS332914"
        sample2.organism = gallus_gallus
        sample2.technology = "RNA-SEQ"
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
        computed_file2.absolute_file_path = "/home/user/data_store/PCL/" + computed_file2.filename
        computed_file2.result = result2
        computed_file2.size_in_bytes = 234
        computed_file2.is_smashable = True
        computed_file2.save()

        assoc2 = SampleComputedFileAssociation()
        assoc2.sample = sample2
        assoc2.computed_file = computed_file2
        assoc2.save()

        dset = Dataset()
        dset.data = {"GSE1487313": ["GSM1487313", "GSM1487222"], "SRX332914": ["SRS332914"]}
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

        self.assertFalse(job.success)

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

            sample = Sample()
            sample.accession_code = file
            sample.title = file
            sample.organism = danio_rerio
            sample.technology = "MICROARRAY"
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
            computed_file.filename = file
            computed_file.absolute_file_path = "/home/user/data_store/raw/TEST/MICROARRAY/" + file
            computed_file.result = result
            computed_file.size_in_bytes = 123
            computed_file.is_smashable = True
            computed_file.save()

            assoc = SampleComputedFileAssociation()
            assoc.sample = sample
            assoc.computed_file = computed_file
            assoc.save()

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

            sample = Sample()
            sample.accession_code = file
            sample.title = file
            sample.organism = danio_rerio
            sample.technology = "RNASEQ"
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
            computed_file.filename = file
            computed_file.absolute_file_path = "/home/user/data_store/raw/TEST/RNASEQ/" + file
            computed_file.result = result
            computed_file.size_in_bytes = 123
            computed_file.is_smashable = True
            computed_file.save()

            assoc = SampleComputedFileAssociation()
            assoc.sample = sample
            assoc.computed_file = computed_file
            assoc.save()

            rnas.append(file)

        # Missing sample that will be filtered
        sample = Sample()
        sample.accession_code = "GSM1487222"
        sample.title = "this sample will be filtered"
        sample.organism = danio_rerio
        sample.technology = "RNASEQ"
        sample.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = experiment
        esa.sample = sample
        esa.save()

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
            final_context["compendium_result"].primary_organism.name, final_context["organism_name"]
        )
        self.assertEqual(final_context["compendium_result"].primary_organism.name, "DANIO_RERIO")
        self.assertEqual(final_context["compendium_result"].organisms.count(), 1)

        # check that sample with no computed file was skipped
        self.assertTrue("GSM1487222" in final_context["filtered_samples"])
        self.assertEqual(
            final_context["filtered_samples"]["GSM1487222"]["experiment_accession_code"], "GSE5678"
        )

        # It's maybe not worth asserting this until we're sure the behavior is correct
        # self.assertEqual(final_context['merged_qn'].shape, (9045, 830))

        zf = zipfile.ZipFile(
            final_context["compendium_result"].result.computedfile_set.first().absolute_file_path
        )
        with zf.open("aggregated_metadata.json") as f:
            metadata = json.load(f)

            # Make sure the data were quantile normalized
            self.assertTrue(metadata.get("quantile_normalized"))
