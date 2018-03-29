from django.contrib import admin
from data_refinery_common.models import ( 
    ProcessorJob,
    DownloaderJob,
    SurveyJob

)
admin.site.register(ProcessorJob)
admin.site.register(DownloaderJob)
admin.site.register(SurveyJob)

