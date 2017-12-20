from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_workers import processors


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
# This should use the get_and_configure_logger.
logger = logging.getLogger(__name__)


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
        valid_processors = list(map(lambda x: x.value, list(ProcessorPipeline)))
        if options["job_id"] is None:
            logger.error("You must specify a job ID.")
            return 1
        elif options["job_name"] not in valid_processors:
            logger.error("You must specify a job name.")
            return 1
        elif options["job_name"] is ProcessorPipeline.AFFY_TO_PCL.value:
            processors.affy_to_pcl(options["job_id"])
            return 0
        elif options["job_name"] is ProcessorPipeline.TRANSCRIPTOME_INDEX:
            processors.build_transcriptome_index(options["job_id"])
            return 0
        elif options["job_name"] is ProcessorPipeline.NO_OP.value:
            processors.no_op_processor(options["job_id"])
            return 0
        else:
            logger.error(("A valid job name was specified for job %s but "
                          "no processor function is known to run it."),
                         str(options["job_id"]))
