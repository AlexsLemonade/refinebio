from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_workers.processors.array_express import affy_to_pcl
from data_refinery_workers.processors.transcriptome_index import build_transcriptome_index
from data_refinery_workers.processors.no_op import no_op_processor


logger = get_and_configure_logger(__name__)


# Test this.
class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--job-name",
            type=str,
            help=("The processor job's name. Must be enumerated in job_lookup."))
        parser.add_argument(
            "--job-id",
            type=int,
            help=("The processor job's ID."))

    def handle(self, *args, **options):
        if options["job_id"] is None:
            logger.error("You must specify a job ID.")
            return 1

        try:
            job_type = ProcessorPipeline[options["job_name"]]
        except KeyError:
            logger.error("You must specify a valid job name.")
            return 1

        if job_type is ProcessorPipeline.AFFY_TO_PCL:
            affy_to_pcl(options["job_id"])
        elif job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX:
            build_transcriptome_index(options["job_id"])
        elif job_type is ProcessorPipeline.NO_OP:
            no_op_processor(options["job_id"])
        else:
            logger.error(("A valid job name was specified for job %s with id %d but "
                          "no processor function is known to run it."),
                         options["job_name"],
                         options["job_id"])
            return 1

        return 0
