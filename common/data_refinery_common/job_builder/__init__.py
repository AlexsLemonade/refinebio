from common.data_refinery_common.models import Sample


class JobBuilder:
    def scope(self):
        """ Returns a queryset with all the samples this class
        knows how to create jobs for """
        return Sample.objects

    def next_job(self, samples):
        """ Given a set of samples return the list of jobs that need to
        be queued for them """
        pass

    def get_downloader_jobs(self, samples):
        pass

    def get_processor_jobs(self, samples):
        pass



class JobFactory:
    def __init__(self, builders = []):
        self.builders = builders

    def get_downloader_jobs(self, samples):
        builder = self._get_builder(samples)
        builder.get_downloader_jobs(samples)

    def get_processor_jobs(self, samples):
        builder = self._get_builder(samples)
        builder.get_processor_jobs(samples)

    def _get_builder(self, samples):
        return self.builders[0]

job_factory = JobFactory()
