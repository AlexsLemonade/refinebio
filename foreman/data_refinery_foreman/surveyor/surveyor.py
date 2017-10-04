import traceback
from django.utils import timezone
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SourceNotSupportedError(BaseException):
    pass


def _get_surveyor_for_source(survey_job: SurveyJob):
    """Factory method for ExternalSourceSurveyors."""
    if survey_job.source_type == "ARRAY_EXPRESS":
        return ArrayExpressSurveyor(survey_job)
    else:
        raise SourceNotSupportedError(
            "Source " + survey_job.source_type + " is not supported.")


def _start_job(survey_job: SurveyJob):
    logger.info("Starting Survey Job #%d for source type: %s.",
                survey_job.id,
                survey_job.source_type)

    survey_job.start_time = timezone.now()
    survey_job.replication_started_at = timezone.now()

    # If the end of the replication range is not already set,
    # set it to the current time.
    if survey_job.replication_ended_at is None:
        survey_job.replication_ended_at = timezone.now()

    survey_job.save()


def _end_job(survey_job: SurveyJob, success=True):
    survey_job.success = success
    survey_job.end_time = timezone.now()
    survey_job.save()


def run_job(survey_job: SurveyJob):
    _start_job(survey_job)

    try:
        surveyor = _get_surveyor_for_source(survey_job)
    except SourceNotSupportedError as e:
        logger.error("Unable to run survey job #%d because: %s",
                     survey_job.id,
                     e)

        _end_job(survey_job, False)
        return survey_job

    try:
        job_success = surveyor.survey()
    except Exception as e:
        logger.error("Exception caught while running job #%d with message: %s",
                     survey_job.id,
                     e)
        logger.error(traceback.format_exc())
        job_success = False

    _end_job(survey_job, job_success)
    return survey_job


def test():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="experiment_accession_code",
                                       value="E-MTAB-3050")
    key_value_pair.save()
    run_job(survey_job)
    return


def survey_experiments(experiments_list_file):
    with open(experiments_list_file, "r") as experiments:
        for experiment in experiments:
            survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
            survey_job.save()
            key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                               key="experiment_accession_code",
                                               value=experiment.rstrip())
            key_value_pair.save()
            run_job(survey_job)
