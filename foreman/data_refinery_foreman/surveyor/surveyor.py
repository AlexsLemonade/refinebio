from django.utils import timezone
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_foreman.surveyor.transcriptome_index import TranscriptomeIndexSurveyor
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class SourceNotSupportedError(BaseException):
    pass


def _get_surveyor_for_source(survey_job: SurveyJob):
    """Factory method for ExternalSourceSurveyors."""
    if survey_job.source_type == "ARRAY_EXPRESS":
        return ArrayExpressSurveyor(survey_job)
    if survey_job.source_type == "SRA":
        return SraSurveyor(survey_job)
    if survey_job.source_type == "TRANSCRIPTOME_INDEX":
        return TranscriptomeIndexSurveyor(survey_job)
    else:
        raise SourceNotSupportedError(
            "Source " + survey_job.source_type + " is not supported.")


def _start_job(survey_job: SurveyJob):
    logger.info("Starting Survey Job for source type: %s.",
                survey_job.source_type,
                survey_job=survey_job.id)

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
        logger.error("Unable to run Survey Job because: %s",
                     e,
                     survey_job=survey_job.id)

        _end_job(survey_job, False)
        return survey_job

    try:
        job_success = surveyor.survey()
    except Exception as e:
        logger.exception("Exception caught while running Survey Job.",
                         survey_job=survey_job.id)
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


def survey_ae_experiment(experiment_accession):
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="experiment_accession_code",
                                       value=experiment_accession)
    key_value_pair.save()
    run_job(survey_job)


def survey_sra_experiments(start_accession, end_accession):
    survey_job = SurveyJob(source_type="SRA")
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="start_accession",
                                       value=start_accession)
    key_value_pair.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="end_accession",
                                       value=end_accession)
    key_value_pair.save()
    run_job(survey_job)
