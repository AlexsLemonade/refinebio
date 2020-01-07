import sys

from django.core.management.base import BaseCommand

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


# Test this.
class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--job-name",
            type=str,
            help=("The processor job's name. Must be enumerated in job_lookup."),
        )
        parser.add_argument("--job-id", type=int, help=("The processor job's ID."))

    def handle(self, *args, **options):
        if options["job_id"] is None:
            logger.error("You must specify a job ID.", job_id=options["job_id"])
            sys.exit(1)

        try:
            job_type = ProcessorPipeline[options["job_name"]]
        except KeyError:
            logger.error(
                "You must specify a valid job name.",
                job_name=options["job_name"],
                job_id=options["job_id"],
            )
            sys.exit(1)

        if job_type is ProcessorPipeline.AFFY_TO_PCL:
            from data_refinery_workers.processors.array_express import affy_to_pcl

            affy_to_pcl(options["job_id"])
        elif job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT:
            from data_refinery_workers.processors.transcriptome_index import (
                build_transcriptome_index,
            )

            build_transcriptome_index(options["job_id"], length="short")
        elif job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG:
            from data_refinery_workers.processors.transcriptome_index import (
                build_transcriptome_index,
            )

            build_transcriptome_index(options["job_id"], length="long")
        elif job_type is ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL:
            from data_refinery_workers.processors.agilent_twocolor import agilent_twocolor_to_pcl

            agilent_twocolor_to_pcl(options["job_id"])
        elif job_type is ProcessorPipeline.ILLUMINA_TO_PCL:
            from data_refinery_workers.processors.illumina import illumina_to_pcl

            illumina_to_pcl(options["job_id"])
        elif job_type is ProcessorPipeline.SALMON:
            from data_refinery_workers.processors.salmon import salmon

            salmon(options["job_id"])
        elif job_type is ProcessorPipeline.TXIMPORT:
            from data_refinery_workers.processors.tximport import tximport

            tximport(options["job_id"])
        elif job_type is ProcessorPipeline.SMASHER:
            from data_refinery_workers.processors.smasher import smash

            smash(options["job_id"])
        elif job_type is ProcessorPipeline.CREATE_COMPENDIA:
            from data_refinery_workers.processors.create_compendia import create_compendia

            create_compendia(options["job_id"])
        elif job_type is ProcessorPipeline.CREATE_QUANTPENDIA:
            from data_refinery_workers.processors.create_quantpendia import create_quantpendia

            create_quantpendia(options["job_id"])
        elif job_type is ProcessorPipeline.NO_OP:
            from data_refinery_workers.processors.no_op import no_op_processor

            no_op_processor(options["job_id"])
        elif job_type is ProcessorPipeline.JANITOR:
            from data_refinery_workers.processors.janitor import run_janitor

            run_janitor(options["job_id"])
        elif job_type is ProcessorPipeline.QN_REFERENCE:
            from data_refinery_workers.processors import qn_reference

            qn_reference.create_qn_reference(options["job_id"])
        else:
            logger.error(
                (
                    "A valid job name was specified for job %s with id %d but "
                    "no processor function is known to run it."
                ),
                options["job_name"],
                options["job_id"],
            )
            sys.exit(1)

        sys.exit(0)
