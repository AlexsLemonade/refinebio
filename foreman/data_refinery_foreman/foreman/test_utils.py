from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)


def create_downloader_job():
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
    og_file_1 = OriginalFile()
    og_file_1.source_filename = "doesn't matter"
    og_file_1.filename = "this either"
    og_file_1.absolute_file_path = "nor this"
    og_file_1.save()

    og_file_2 = OriginalFile()
    og_file_2.source_filename = "doesn't matter"
    og_file_2.filename = "this either"
    og_file_2.absolute_file_path = "nor this"
    og_file_2.save()

    if pipeline == "AFFY_TO_PCL":
        downloader_job = DownloaderJob(
            downloader_task="SRA",
            batch_job_id="DEFAULT",
            num_retries=0,
            accession_code="NUNYA",
            success=None,
        )
        downloader_job.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = og_file_2
        assoc.downloader_job = downloader_job
        assoc.save()

        assoc1 = DownloaderJobOriginalFileAssociation()
        assoc1.original_file = og_file_1
        assoc1.downloader_job = downloader_job
        assoc1.save()

    processor_job = ProcessorJob(
        downloader_job=downloader_job,
        pipeline_applied=pipeline,
        batch_job_id="PROCESSOR/dispatch-1528945054-e8eaf540",
        ram_amount=ram_amount,
        num_retries=0,
        success=None,
        start_time=start_time,
    )
    processor_job.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file_1
    assoc1.processor_job = processor_job
    assoc1.save()

    assoc = ProcessorJobOriginalFileAssociation()
    assoc.original_file = og_file_2
    assoc.processor_job = processor_job
    assoc.save()

    return processor_job


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
