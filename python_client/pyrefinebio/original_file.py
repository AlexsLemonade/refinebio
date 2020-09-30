import pyrefinebio.job as prb_job
import pyrefinebio.sample as prb_sample
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class OriginalFile:
    def __init__(
        self,
        id=None,
        filename=None,
        size_in_bytes=None,
        sha1=None,
        samples=[],
        processor_jobs=[],
        downloader_jobs=[],
        source_url=None,
        source_filename=None,
        is_downloaded=None,
        is_archive=None,
        has_raw=None,
        created_at=None,
        last_modified=None,
    ):
        self.id = id
        self.filename = filename
        self.size_in_bytes = size_in_bytes
        self.sha1 = sha1
        self.samples = [prb_sample.Sample(**sample) for sample in samples]
        self.processor_jobs = [
            prb_job.ProcessorJob(**processor_job) for processor_job in processor_jobs
        ]
        self.downloader_jobs = [
            prb_job.DownloaderJob(**downloader_job) for downloader_job in downloader_jobs
        ]
        self.source_url = source_url
        self.source_filename = source_filename
        self.is_downloaded = is_downloaded
        self.is_archive = is_archive
        self.has_raw = has_raw
        self.created_at = created_at
        self.last_modified = last_modified

    @classmethod
    def get(cls, id):
        response = get_by_endpoint("original_files/" + str(id))
        return OriginalFile(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("original_files", params=kwargs)
        return generator_from_pagination(response, cls)
