import datetime
import traceback
from data_refinery_models.models import SurveyJob, SurveyJobKeyValue
from .array_express_surveyor import ArrayExpressSurveyor


def SourceNotSupportedError(Exception):
    pass


def _get_surveyor_for_source(survey_job: SurveyJob):
    """Factory method for ExternalSourceSurveyors."""
    if(survey_job.source_type == "ARRAY_EXPRESS"):
        return ArrayExpressSurveyor(survey_job)
    else:
        raise SourceNotSupportedError(
            "Source " + survey_job.source_type + " is not supported.")


def _start_job(survey_job: SurveyJob):
    survey_job.start_time = datetime.datetime.now(datetime.timezone.utc)

    if survey_job.replication_ended_at is None:
        survey_job.replication_ended_at = datetime.datetime.now(
            datetime.timezone.utc)

    survey_job.save()


def _end_job(survey_job: SurveyJob, success=True):
    survey_job.success = success
    survey_job.end_time = datetime.datetime.now(datetime.timezone.utc)
    survey_job.save()


def run_job(survey_job: SurveyJob):
    _start_job(survey_job)

    try:
        surveyor = _get_surveyor_for_source(survey_job)
    except SourceNotSupportedError as e:
        # This should be logging, not printing. I need to set that up.
        log_message = "Unable to run survey job # " + survey_job.id
        log_message = log_message + " because: " + str(e)
        print(log_message)

        _end_job(survey_job, False)
        return survey_job

    try:
        job_success = surveyor.survey(survey_job)
    except Exception as e:
        print("Exception caught while running job #" + str(survey_job.id) +
              " with message: " + str(e))
        print(traceback.format_exc())
        job_success = False

    _end_job(survey_job, job_success)
    return survey_job


def test():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="accession_code",
                                       value="A-AFFY-1")
    key_value_pair.save()
    run_job(survey_job)
    return
