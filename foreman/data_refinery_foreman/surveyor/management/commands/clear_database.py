"""
This command will create and run survey jobs for the specified ensembl division.
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import *

class Command(BaseCommand):
    def handle(self, *args, **options):
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()
        Experiment.objects.all().delete()
        ExperimentAnnotation.objects.all().delete()
        ComputationalResult.objects.all().delete()
        ComputationalResultAnnotation.objects.all().delete()
        OrganismIndex.objects.all().delete()
        OriginalFile.objects.all().delete()
        ComputedFile.objects.all().delete()
        ExperimentSampleAssociation.objects.all().delete()
        DownloaderJobOriginalFileAssociation.objects.all().delete()
        ProcessorJobOriginalFileAssociation.objects.all().delete()
        SampleResultAssociation.objects.all().delete()
        Organism.objects.all().delete()
        SurveyJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        DownloaderJob.objects.all().delete()
        ProcessorJob.objects.all().delete()
