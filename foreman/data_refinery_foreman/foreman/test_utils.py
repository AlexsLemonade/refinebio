from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)


def create_downloader_job(suffix="e8eaf540"):
    job = DownloaderJob(
        downloader_task="SRA",
        batch_job_id="DEFAULT",
        num_retries=0,
        accession_code="NUNYA",
        success=None,
    )
    job.save()

    og_file = OriginalFile()
    og_file.source_filename = "doesn't matter"
    og_file.filename = "this either"
    og_file.absolute_file_path = "nor this"
    og_file.save()

    assoc1 = DownloaderJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.downloader_job = job
    assoc1.save()

    og_file = OriginalFile()
    og_file.source_filename = "doesn't matter"
    og_file.filename = "this either"
    og_file.absolute_file_path = "nor this"
    og_file.save()

    assoc = DownloaderJobOriginalFileAssociation()
    assoc.original_file = og_file
    assoc.downloader_job = job
    assoc.save()

    return job


def create_processor_job(pipeline="AFFY_TO_PCL", ram_amount=2048, start_time=None):
    job = ProcessorJob(
        pipeline_applied=pipeline,
        batch_job_id="PROCESSOR/dispatch-1528945054-e8eaf540",
        ram_amount=ram_amount,
        num_retries=0,
        volume_index="1",
        success=None,
        start_time=start_time,
    )
    job.save()

    og_file = OriginalFile()
    og_file.source_filename = "doesn't matter"
    og_file.filename = "this either"
    og_file.absolute_file_path = "nor this"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = job
    assoc1.save()

    og_file = OriginalFile()
    og_file.source_filename = "doesn't matter"
    og_file.filename = "this either"
    og_file.absolute_file_path = "nor this"
    og_file.save()

    assoc = ProcessorJobOriginalFileAssociation()
    assoc.original_file = og_file
    assoc.processor_job = job
    assoc.save()

    return job


def create_survey_job():
    job = SurveyJob(
        source_type="SRA",
        batch_job_id="SURVEYOR/dispatch-1528945054-e8eaf540",
        num_retries=0,
        success=None,
    )

    job.save()

    sjkv = SurveyJobKeyValue()
    sjkv.key = "experiment_accession_code"
    sjkv.value = "RJ-1234-XYZ"
    sjkv.survey_job = job
    sjkv.save()

    return job
