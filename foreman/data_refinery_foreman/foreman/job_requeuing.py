from data_refinery_common.job_lookup import (
    SMASHER_JOB_TYPES,
    Downloaders,
    ProcessorPipeline,
    SurveyJobTypes,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)

logger = get_and_configure_logger(__name__)


def requeue_downloader_job(last_job: DownloaderJob) -> (bool, str):
    """Queues a new downloader job.

    The new downloader job will have num_retries one greater than
    last_job.num_retries.

    Returns True and the volume index of the downloader job upon successful dispatching,
    False and an empty string otherwise.
    """
    num_retries = last_job.num_retries + 1

    ram_amount = last_job.ram_amount
    # If there's no start time then it's likely that the instance got
    # cycled which means we didn't get OOM-killed, so we don't need to
    # increase the RAM amount.
    if last_job.start_time and last_job.failure_reason is None:
        if ram_amount == 1024:
            ram_amount = 4096
        elif ram_amount == 4096:
            ram_amount = 16384

    original_file = last_job.original_files.first()

    if not original_file:
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = (
            "Foreman told to requeue a DownloaderJob without an OriginalFile - why?!"
        )
        last_job.save()
        logger.info(
            "Foreman told to requeue a DownloaderJob without an OriginalFile - why?!",
            last_job=str(last_job),
        )
        return False, ""

    if not original_file.needs_processing():
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = "Foreman told to redownload job with prior successful processing."
        last_job.save()
        logger.info(
            "Foreman told to redownload job with prior successful processing.",
            last_job=str(last_job),
        )
        return False, ""

    first_sample = original_file.samples.first()

    # This is a magic string that all the dbGaP studies appear to have
    if first_sample and ("in the dbGaP study" in first_sample.title):
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = "Sample is dbGaP access controlled."
        last_job.save()
        logger.info(
            "Avoiding requeuing for DownloaderJob for dbGaP run accession: "
            + str(first_sample.accession_code)
        )
        return False, ""

    new_job = DownloaderJob(
        num_retries=num_retries,
        downloader_task=last_job.downloader_task,
        ram_amount=ram_amount,
        accession_code=last_job.accession_code,
        was_recreated=last_job.was_recreated,
        # For now default to 0, will need to be calculated.
        volume_index=0,
    )
    new_job.save()

    for original_file in last_job.original_files.all():
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=new_job, original_file=original_file
        )

    logger.debug(
        "Requeuing Downloader Job which had ID %d with a new Downloader Job with ID %d.",
        last_job.id,
        new_job.id,
    )
    try:
        if send_job(Downloaders[last_job.downloader_task], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            new_job.delete()
            return False, ""
    except Exception:
        logger.error(
            "Failed to requeue DownloaderJob which had ID %d with a new DownloaderJob with ID %d.",
            last_job.id,
            new_job.id,
        )
        # Can't communicate with Batch just now, leave the job for a later loop.
        new_job.delete()
        return False, ""

    return True, new_job.volume_index


def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    # The Salmon pipeline is quite RAM-sensitive.
    # Try it again with an increased RAM amount, if possible.
    new_ram_amount = last_job.ram_amount

    # If there's no start time then it's likely that the instance got
    # cycled which means we didn't get OOM-killed, so we don't need to
    # increase the RAM amount.
    if last_job.start_time:
        # There's only one size of tximport jobs.
        if last_job.pipeline_applied == "TXIMPORT":
            new_ram_amount = 32768
        # These initial values are set in common/job_lookup.py:determine_ram_amount
        elif last_job.pipeline_applied == "SALMON" or last_job.pipeline_applied.startswith(
            "TRANSCRIPTOME"
        ):
            if new_ram_amount == 4096:
                new_ram_amount = 8192
            if new_ram_amount == 8192:
                new_ram_amount = 12288
            elif new_ram_amount == 12288:
                new_ram_amount = 16384
            elif new_ram_amount == 16384:
                new_ram_amount = 32768
            elif new_ram_amount == 32768:
                new_ram_amount = 65536
        # The AFFY pipeline is somewhat RAM-sensitive.
        # Also NO_OP can fail and be retried, so we want to attempt ramping up ram.
        # Try it again with an increased RAM amount, if possible.
        elif last_job.pipeline_applied == "AFFY_TO_PCL" or last_job.pipeline_applied == "NO_OP":
            if new_ram_amount == 2048:
                new_ram_amount = 4096
            elif new_ram_amount == 4096:
                new_ram_amount = 8192
            elif new_ram_amount == 8192:
                new_ram_amount = 32768
        elif (
            last_job.pipeline_applied == "ILLUMINA_TO_PCL"
            and "non-zero exit status -9" in last_job.failure_reason
        ):
            if new_ram_amount == 2048:
                new_ram_amount = 4096
            elif new_ram_amount == 4096:
                new_ram_amount = 8192

    volume_index = last_job.volume_index
    # Make sure volume_index is set to something, unless it's a
    # smasher job type because the smasher instance doesn't have a
    # volume_index.
    if (not volume_index or volume_index == "-1") and ProcessorPipeline[
        last_job.pipeline_applied
    ] not in SMASHER_JOB_TYPES:
        # Default to 0 for now. At some point we'll have to detect
        # when 0 is full, and if so go to 1, etc.
        volume_index = 0

    new_job = ProcessorJob(
        num_retries=num_retries,
        pipeline_applied=last_job.pipeline_applied,
        ram_amount=new_ram_amount,
        volume_index=volume_index,
    )
    new_job.save()

    for original_file in last_job.original_files.all():
        ProcessorJobOriginalFileAssociation.objects.get_or_create(
            processor_job=new_job, original_file=original_file
        )

    for dataset in last_job.datasets.all():
        ProcessorJobDatasetAssociation.objects.get_or_create(processor_job=new_job, dataset=dataset)

    try:
        logger.debug(
            "Requeuing Processor Job which had ID %d with a new Processor Job with ID %d.",
            last_job.id,
            new_job.id,
        )
        if send_job(ProcessorPipeline[last_job.pipeline_applied], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            new_job.delete()
    except Exception:
        logger.warn(
            "Failed to requeue Processor Job which had ID %d with a new Processor Job with ID %d.",
            last_job.id,
            new_job.id,
            exc_info=1,
        )
        # Can't communicate with Batch just now, leave the job for a later loop.
        new_job.delete()


def requeue_survey_job(last_job: SurveyJob) -> None:
    """Queues a new survey job.

    The new survey job will have num_retries one greater than
    last_job.num_retries.
    """

    num_retries = last_job.num_retries + 1

    new_job = SurveyJob(num_retries=num_retries, source_type=last_job.source_type)

    if new_job.num_retries == 1:
        new_job.ram_amount = 4096
    elif new_job.num_retries in [2, 3]:
        new_job.ram_amount = 16384
    else:
        new_job.ram_amount = 1024

    new_job.save()

    keyvalues = SurveyJobKeyValue.objects.filter(survey_job=last_job)

    for keyvalue in keyvalues:
        SurveyJobKeyValue.objects.get_or_create(
            survey_job=new_job, key=keyvalue.key, value=keyvalue.value,
        )

    logger.debug(
        "Requeuing SurveyJob which had ID %d with a new SurveyJob with ID %d.",
        last_job.id,
        new_job.id,
    )

    try:
        if send_job(SurveyJobTypes.SURVEYOR, job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            new_job.delete()
    except Exception:
        logger.error(
            "Failed to requeue Survey Job which had ID %d with a new Surevey Job with ID %d.",
            last_job.id,
            new_job.id,
        )
        # Can't communicate with AWS just now, leave the job for a later loop.
        new_job.delete()

    return True
