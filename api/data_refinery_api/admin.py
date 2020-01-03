from django.contrib import admin
from data_refinery_common.models import (
    ProcessorJob,
    DownloaderJob,
    SurveyJob,
    Experiment,
    Sample,
    OriginalFile,
    ComputationalResult,
    OrganismIndex,
)

admin.site.register(ProcessorJob)
admin.site.register(DownloaderJob)
admin.site.register(SurveyJob)

admin.site.register(Experiment)
admin.site.register(Sample)
admin.site.register(OriginalFile)
admin.site.register(ComputationalResult)
admin.site.register(OrganismIndex)
