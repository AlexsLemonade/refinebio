"""This command is intended for development purposes.
It creates the database records necessary for a processor job to
run and queues one. It assumes that the files
/home/user/data_store/raw/IlluminaHiSeq2500/SALMON/ERR1680082_1.fastq
and
/home/user/data_store/raw/IlluminaHiSeq2500/SALMON/ERR1680082_2.fastq
exist. The easiest way to run this is with the tester.sh script.
(Changing queue_downloader to queue_processor.)
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchKeyValue,
    File,
    BatchStatuses,
    ProcessorJob
)

from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("processor-type",
                            help=("The type of processor. Valid options are:\n"
                                  "[SRA, TRANSCRIPTOME_INDEX]"))

    def run_sra_processor(self):
        # Create all the dummy data that would have been created
        # before a processor job could have been generated.
        survey_job = SurveyJob(
            source_type="SRA"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="SRA",
            pipeline_required="SALMON",
            platform_accession_code="IlluminaHiSeq2500",
            experiment_accession_code="PRJEB5018",
            experiment_title="It doesn't really matter.",
            organism_id=10090,
            organism_name="MUS MUSCULUS",
            release_date="2014-03-25",
            last_uploaded_date="2016-05-20",
            status=BatchStatuses.NEW.value,
        )
        batch.save()

        File(name="ERR1680082_1.fastq",
             download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/"
                           "ERR1680082/ERR1680082_1.fastq.gz"),
             raw_format="fastq",
             processed_format="sf",
             internal_location="IlluminaHiSeq2500/SALMON",
             size_in_bytes=2214725074,
             batch=batch).save()

        File(name="ERR1680082_2.fastq",
             download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/"
                           "ERR1680082/ERR1680082_2.fastq.gz"),
             raw_format="fastq",
             processed_format="sf",
             internal_location="IlluminaHiSeq2500/SALMON",
             size_in_bytes=2214725074,
             batch=batch).save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        logger.info("Queuing a processor job.")
        send_job(ProcessorPipeline[batch.pipeline_required], processor_job.id)

    def run_trasnscriptome_processor(self):
        # Create all the dummy data that would have been created
        # before a processor job could have been generated.
        survey_job = SurveyJob(
            source_type="TRANSCRIPTOME_INDEX"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="TRANSCRIPTOME_INDEX",
            pipeline_required="TRANSCRIPTOME_INDEX",
            platform_accession_code="EnsemblPlants",
            experiment_accession_code="aegilops_tauschii",
            experiment_title="It doesn't really matter.",
            organism_id=37682,
            organism_name="AEGILOPS TAUSCHII",
            release_date="2017-11-02",
            last_uploaded_date="2017-11-02",
            status=BatchStatuses.DOWNLOADED.value,
        )
        batch.save()

        kmer_size_property = BatchKeyValue(batch=batch,
                                           key="kmer_size",
                                           value="31")
        kmer_size_property.save()

        gtf_file = File(name="aegilops_tauschii_short.gtf.gz",
                        download_url=("ftp://ftp.ensemblgenomes.org/pub/release-37/plants/gtf"
                                      "/aegilops_tauschii/Aegilops_tauschii.ASM34733v1.37.gtf.gz"),
                        raw_format="gtf.gz",
                        processed_format="tar.gz",
                        internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
                        size_in_bytes=-1,
                        batch=batch)
        gtf_file.save()

        fasta_file = File(name="aegilops_tauschii_short.fa.gz",
                          download_url=("ftp://ftp.ensemblgenomes.org/pub/release-37/plants/fasta"
                                        "/aegilops_tauschii/dna/Aegilops_tauschii."
                                        "ASM34733v1.dna.toplevel.fa.gz"),
                          raw_format="fa.gz",
                          processed_format="tar.gz",
                          internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
                          size_in_bytes=-1,
                          batch=batch)
        fasta_file.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        logger.info("Queuing a processor job.")
        send_job(ProcessorPipeline[batch.pipeline_required], processor_job.id)

    def handle(self, *args, **options):
        if options["processor-type"] is None:
            logger.error("You must set the type of processor to queue.")
            return 1
        elif options["processor-type"] == "SRA":
            self.run_sra_processor()
            return 0
        elif options["processor-type"] == "TRANSCRIPTOME_INDEX":
            self.run_trasnscriptome_processor()
            return 0
        else:
            logger.error("Unrecognized processor-type.")
            return 1
