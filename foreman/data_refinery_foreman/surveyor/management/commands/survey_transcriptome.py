"""
This command will create and run survey jobs for the specified ensembl division.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "ensembl_division",
            help=("An ensembl division. Possible values are:"
                  "['Ensembl', 'EnsemblFungi', 'EnsemblBacteria', 'EnsemblProtists',"
                  "'EnsemblMetazoa', 'EnsemblPlants']"))
        parser.add_argument(
            "number_of_organisms",
            nargs='?',
            default="-1",
            help=("How many organisms to survey in this run. Omitting this argument or passing"
                  " -1 will survey all organisms in the division."))

    def handle(self, *args, **options):
        if options["ensembl_division"] is None:
            logger.error("You must specify an ensembl_division.")
            return 1
        else:
            survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
            survey_job.save()
            key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                               key="ensembl_division",
                                               value=options["ensembl_division"])
            key_value_pair.save()

            key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                               key="number_of_organisms",
                                               value=options["number_of_organisms"])
            key_value_pair.save()

            surveyor.run_job(survey_job)
            return 0
