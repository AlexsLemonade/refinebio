import signal
import sys

from django.utils import timezone

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor
from data_refinery_foreman.surveyor.geo import GeoSurveyor
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_foreman.surveyor.transcriptome_index import TranscriptomeIndexSurveyor


logger = get_and_configure_logger(__name__)

CURRENT_JOB = None

class SourceNotSupportedError(Exception):
    pass

def signal_handler(sig, frame):
    """Signal Handler, works for both SIGTERM and SIGINT"""
    global CURRENT_JOB
    if CURRENT_JOB:
        CURRENT_JOB.success = False
        CURRENT_JOB.end_time = timezone.now()
        CURRENT_JOB.num_retries = CURRENT_JOB.num_retries - 1
        CURRENT_JOB.failure_reason = "Interruped by SIGTERM/SIGINT: " + str(sig)
        CURRENT_JOB.save()

    sys.exit(0)

def _get_surveyor_for_source(survey_job: SurveyJob):
    """Factory method for ExternalSourceSurveyors."""
    if survey_job.source_type == "ARRAY_EXPRESS":
        return ArrayExpressSurveyor(survey_job)
    if survey_job.source_type == "SRA":
        return SraSurveyor(survey_job)
    if survey_job.source_type == "TRANSCRIPTOME_INDEX":
        return TranscriptomeIndexSurveyor(survey_job)
    if survey_job.source_type == "GEO":
        return GeoSurveyor(survey_job)
    else:
        raise SourceNotSupportedError(
            "Source " + survey_job.source_type + " is not supported.")


def _start_job(survey_job: SurveyJob) -> SurveyJob:
    """Start survey job, setting time properties."""
    logger.debug("Starting Survey Job for source type: %s.",
                survey_job.source_type,
                survey_job=survey_job.id)

    # Set up the SIGTERM handler so we can appropriately handle being interrupted.
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    survey_job.start_time = timezone.now()
    survey_job.save()

    global CURRENT_JOB
    CURRENT_JOB = survey_job

    return survey_job


def _end_job(survey_job: SurveyJob, success=True) -> SurveyJob:
    """Ends survey job, setting success and time properties."""
    survey_job.success = success
    survey_job.end_time = timezone.now()
    survey_job.save()

    return survey_job


def run_job(survey_job: SurveyJob) -> SurveyJob:
    """Runs a survey job and handles errors."""
    survey_job = _start_job(survey_job)

    try:
        surveyor = _get_surveyor_for_source(survey_job)
    except SourceNotSupportedError as e:
        logger.error("Unable to run Survey Job because: %s",
                     e,
                     survey_job=survey_job.id)

        _end_job(survey_job, False)
        return survey_job

    try:
        job_success = surveyor.survey(source_type=survey_job.source_type)
    except Exception as e:
        logger.exception("Exception caught while running Survey Job.",
                         survey_job=survey_job.id)
        survey_job.failure_reason = str(e)
        job_success = False

    survey_job = _end_job(survey_job, job_success)
    return survey_job


def survey_experiment(experiment_accession: str, source_type: str):
    """Survey an experiment of type `source_type`.

    Source type corresponds to one of the external sources we
    support. It must be one of the following values:
      * SRA
      * GEO
      * ARRAY_EXPRESS
    """
    survey_job = SurveyJob(source_type=source_type)
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="experiment_accession_code",
                                       value=experiment_accession)
    key_value_pair.save()
    run_job(survey_job)

    return survey_job


def survey_transcriptome_index(organism_name=None, ensembl_division='Ensembl'):
    """Special one-off surveyor to build transcriptome indices.

    The external source this uses is ensembl.org which is divided into
    multiple divisions. This function surveys only one division at a
    time. If an `organism_name` is provided, survey only that
    organism, otherwise survey the entire division.
    """
    survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
    survey_job.save()
    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="ensembl_division",
                                       value=ensembl_division)
    key_value_pair.save()

    if organism_name:
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="organism_name",
                                           value=organism_name)
        key_value_pair.save()

    run_job(survey_job)

    return survey_job
