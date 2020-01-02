from django.test import TransactionTestCase
from django.core.management import call_command

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    ExperimentOrganismAssociation,
)


class QuantendiaCommandTestCase(TransactionTestCase):
    def test_quantpendia_command(self):
        organism = get_organism_with_qn_target()
        make_test_data(organism)

        try:
            call_command("create_quantpendia", organisms=organism.name)
        except SystemExit as e:  # this is okay!
            pass

        processor_job = (
            ProcessorJob.objects.filter(pipeline_applied="CREATE_QUANTPENDIA")
            .order_by("-created_at")
            .first()
        )

        # check that the processor job was created correctly
        self.assertIsNotNone(processor_job)
        self.assertEquals(processor_job.datasets.first().data, {"GSE51088": ["GSM1237818"]})


def get_organism_with_qn_target():
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
    return danio_rerio


def make_test_data(organism):
    experiment = Experiment()
    experiment.accession_code = "GSE51088"
    experiment.technology = "RNA-SEQ"
    experiment.save()

    xoa = ExperimentOrganismAssociation()
    xoa.experiment = experiment
    xoa.organism = organism
    xoa.save()

    result = ComputationalResult()
    result.save()

    sample = Sample()
    sample.accession_code = "GSM1237818"
    sample.title = "GSM1237818"
    sample.organism = organism
    sample.technology = "RNA-SEQ"
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
