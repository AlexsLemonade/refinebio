from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class DownloaderJob:
    def __init__(
        self,
        id=None,
        downloader_task=None,
        num_retries=None,
        retried=None,
        was_recorded=None,
        worker_id=None,
        worker_version=None,
        nomad_job_id=None,
        failure_reason=None,
        success=None,
        original_files=None,
        start_time=None,
        end_time=None,
        created_at=None,
        last_modified=None,
    ):
        self.id = id
        self.downloader_task = downloader_task
        self.num_retries = num_retries
        self.retried = retried
        self.was_recorded = was_recorded
        self.worker_id = worker_id
        self.worker_version = worker_version
        self.nomad_job_id = nomad_job_id
        self.failure_reason = failure_reason
        self.success = success
        self.original_files = original_files
        self.start_time = start_time
        self.end_time = end_time
        self.created_at = created_at
        self.last_modified = last_modified

    @classmethod
    def get(cls, id):
        response = get_by_endpoint("jobs/downloader/" + str(id))
        return DownloaderJob(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("jobs/downloader", params=kwargs)
        return generator_from_pagination(response, cls)


class ProcessorJob:
    def __init__(
        self,
        id=None,
        pipeline_applied=None,
        num_retries=None,
        retried=None,
        worker_id=None,
        ram_amount=None,
        volume_index=None,
        worker_version=None,
        failure_reason=None,
        nomad_job_id=None,
        success=None,
        original_files=None,
        datasets=None,
        start_time=None,
        end_time=None,
        created_at=None,
        last_modified=None,
    ):
        self.id = id
        self.pipeline_applied = pipeline_applied
        self.num_retries = num_retries
        self.retried = retried
        self.worker_id = worker_id
        self.ram_amount = ram_amount
        self.volume_index = volume_index
        self.worker_version = worker_version
        self.failure_reason = failure_reason
        self.nomad_job_id = nomad_job_id
        self.success = success
        self.original_files = original_files
        self.datasets = datasets
        self.start_time = start_time
        self.end_time = end_time
        self.created_at = created_at
        self.last_modified = last_modified

    @classmethod
    def get(cls, id):
        response = get_by_endpoint("jobs/processor/" + str(id))
        return ProcessorJob(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("jobs/processor", params=kwargs)
        return generator_from_pagination(response, cls)


class SurveyJob:
    def __init__(
        self,
        id=None,
        source_type=None,
        success=None,
        start_time=None,
        end_time=None,
        created_at=None,
        last_modified=None,
    ):
        self.id = id
        self.source_type = source_type
        self.success = success
        self.start_time = start_time
        self.end_time = end_time
        self.created_at = created_at
        self.last_modified = last_modified

    @classmethod
    def get(cls, id):
        response = get_by_endpoint("jobs/survey/" + str(id))
        return SurveyJob(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("jobs/survey", params=kwargs)
        return generator_from_pagination(response, cls)
