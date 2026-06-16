"""Calculate RAM-hours consumed processing a sample or experiment.

A RAM-hour is 1 GB of RAM held for 1 hour. Useful for understanding how
expensive a given sample/experiment was to process, especially when deciding
whether to manually debug a stuck or repeated job vs. let it ride.

All downloader and processor jobs for the sample are counted regardless of
their success/failure status. Only jobs with both start_time and end_time set
contribute (in-flight jobs are skipped — they haven't finished yet).

Invocation:
  Local dev:  rbio compose:manage foreman calculate_ram_hours GSE12345
  Prod:       ssh <foreman-ec2>; docker exec dr_foreman \\
                  python3 manage.py calculate_ram_hours GSE12345
"""

import csv
import sys
from datetime import datetime, timezone

from django.core.management.base import BaseCommand

from data_refinery_common.models import Experiment, Sample


class Command(BaseCommand):
    help = __doc__

    def add_arguments(self, parser):
        parser.add_argument(
            "accession_code",
            help=(
                "Experiment accession by default, or sample accession with --sample. "
                "For transcriptome indices, an organism name like 'Homo Sapiens'."
            ),
        )
        parser.add_argument(
            "--start-date",
            default="2021-06-23",
            type=lambda s: datetime.strptime(s, "%Y-%m-%d").replace(tzinfo=timezone.utc),
            help=(
                "Skip jobs that started before this date (YYYY-MM-DD). "
                "Useful for scoping to a specific reprocessing run. Default: 2021-06-23."
            ),
        )
        parser.add_argument(
            "--sample",
            action="store_true",
            help=(
                "Treat accession_code as a sample accession rather than an experiment. "
                "Otherwise the accession is resolved as an experiment and every sample "
                "in it is included."
            ),
        )

    def handle(self, *args, **options):
        accession = options["accession_code"]
        start_date = options["start_date"]
        is_sample = options["sample"]

        writer = csv.writer(sys.stdout)
        writer.writerow(
            ["accession_code", "downloader_ram_hours", "processor_ram_hours", "total_ram_hours"]
        )

        if is_sample:
            self._emit_sample(writer, accession, start_date)
        else:
            try:
                experiment = Experiment.objects.get(accession_code=accession)
            except Experiment.DoesNotExist:
                self.stderr.write(f"Experiment {accession!r} not found.")
                sys.exit(1)
            for sample in experiment.samples.all():
                self._emit_sample(writer, sample.accession_code, start_date)

    def _emit_sample(self, writer, accession_code, start_date):
        try:
            sample = Sample.objects.get(accession_code=accession_code)
        except Sample.DoesNotExist:
            self.stderr.write(f"Sample {accession_code!r} not found.")
            return

        downloader_ram_hours = sum(
            _ram_hours_for(job)
            for job in sample.get_downloader_jobs()
            if job.start_time and job.start_time > start_date
        )
        processor_ram_hours = sum(
            _ram_hours_for(job)
            for job in sample.get_processor_jobs()
            if job.start_time and job.start_time > start_date
        )
        total = downloader_ram_hours + processor_ram_hours

        writer.writerow([accession_code, downloader_ram_hours, processor_ram_hours, total])


def _ram_hours_for(job):
    """Return RAM-hours (GB * hours) for one job, or 0 if it never ran to completion."""
    if not (job.start_time and job.end_time and job.ram_amount):
        return 0
    hours = (job.end_time - job.start_time).total_seconds() / 3600
    gigabytes = job.ram_amount / 1024  # ram_amount is in MB
    return hours * gigabytes
