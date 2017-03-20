import datetime
from bioinformatics_mill_models.models import SurveyJob
from .array_express_surveyor import ArrayExpressSurveyor


def get_surveyor_for_source(survey_job: SurveyJob):
    """Factory method for ExternalSourceSurveyors."""
    if(survey_job.source_type == "ARRAY_EXPRESS"):
        return ArrayExpressSurveyor(survey_job)
    else:
        # Probably should be more specific, but I'm not sure what exception to
        # use yet
        raise Exception("Source " + source_type + " is not supported.")


def start_job(survey_job: SurveyJob):
    survey_job.start_time = datetime.datetime.now(datetime.timezone.utc)

    if(survey_job.replication_ended_at == None):
        survey_job.replication_ended_at = datetime.datetime.now(
            datetime.timezone.utc)

    survey_job.save()


def end_job(survey_job: SurveyJob, success=True):
    survey_job.success = success
    survey_job.end_time = datetime.datetime.now(datetime.timezone.utc)
    survey_job.save()


def run_job(survey_job: SurveyJob):
    start_job(survey_job)

    try:
        surveyor = get_surveyor_for_source(survey_job)
    # once again, should be more specific
    except Exception as e:
        # This should be logging, not printing. I need to set that up.
        log_message = "Unable to run survey job # " + survey_job.id
        log_message = log_message + " because: " + str(e)
        print(log_message)

        end_job(survey_job, False)
        return survey_job

    job_success = surveyor.survey(survey_job)
    end_job(survey_job, job_success)
    return survey_job


def go():
    run_job(SurveyJob(source_type="ARRAY_EXPRESS"))
    return
