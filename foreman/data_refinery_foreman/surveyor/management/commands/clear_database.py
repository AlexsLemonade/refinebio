"""
This command will clear out the database to make repeating tests easier.
"""

from django.core.management.base import BaseCommand

from data_refinery_common.models import *


class Command(BaseCommand):
    def handle(self, *args, **options):
        print("-------------------------------------------------------------------------------")
        print("This will delete all objects in the database. Are you sure you want to do this?")
        answer = input('You must type "yes", all other input will be ignored: ')

        if answer == "yes":
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

            print("The database has been cleared.")
        else:
            print("The database has been left alone.")
