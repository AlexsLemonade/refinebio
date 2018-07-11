"""
This command will clear out the database to make repeating tests easier.
"""

from django.core.management.base import BaseCommand
from datetime import datetime
import subprocess
import os
import json
from data_refinery_common.models import OriginalFile, ProcessorJob,
ProcessorJobOriginalFileAssociation, Sample, OriginalFileSampleAssociation

external_data_directory = os.curdir + "/data_store/external-data/"


class Command(BaseCommand):
    def handle(self, *args, **options):
        affy_jobs = []

        for filename in os.listdir(external_data_directory):
            if filename.endswith(".CEL"):
                job = ProcessorJob()
                job.pipeline_applied = "AFFY_TO_PCL"
                job.save()

                file = OriginalFile()
                file.source_filename = filename
                file.filename = filename
                file.absolute_file_path = os.path.abspath(external_data_directory + filename)
                file.save()

                og_file_assoc = ProcessorJobOriginalFileAssociation()
                og_file_assoc.original_file = file
                og_file_assoc.processor_job = job
                og_file_assoc.save()

                sample = Sample()
                sample.accession_code = str(datetime.now().timestamp())
                sample.save()

                sample_assoc = OriginalFileSampleAssociation()
                sample_assoc.sample = sample
                sample_assoc.original_file = file
                sample_assoc.save()

                affy_jobs.append(job.id)

                continue
            else:
                continue
        with open('data_store/processor-jobs.txt', 'w') as file:
            file.write(json.dumps({"affy": affy_jobs}))
