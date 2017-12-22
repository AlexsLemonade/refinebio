from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_workers.downloaders.array_express import download_array_express
from data_refinery_workers.downloaders.transcriptome_index import download_transcriptome
from data_refinery_workers.downloaders.sra import download_sra


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--job-name",
            type=str,
            help=("The downloader job's name. Must be enumerated in job_lookup."))
        parser.add_argument(
            "--job-id",
            type=int,
            help=("The downloader job's ID."))

    def handle(self, *args, **options):
        if options["job_id"] is None:
            logger.error("You must specify a job ID.")
            return 1

        try:
            job_type = Downloaders[options["job_name"]]
        except KeyError:
            logger.error("You must specify a valid job name.")
            return 1

        if job_type is Downloaders.ARRAY_EXPRESS:
            download_array_express(options["job_id"])
        elif job_type is Downloaders.TRANSCRIPTOME_INDEX:
            download_transcriptome(options["job_id"])
        elif job_type is Downloaders.SRA:
            download_sra(options["job_id"])
        else:
            logger.error(("A valid job name was specified for job %s with id %d but "
                          "no downloader function is known to run it."),
                         options["job_name"],
                         options["job_id"])
            return 1

        return 0
