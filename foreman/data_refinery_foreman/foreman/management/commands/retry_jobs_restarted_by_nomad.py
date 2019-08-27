# Find the samples where the latest job that was applied to them
# failed with because Nomad restarted it, and requeue downloader jobs
# for them
# ref https://github.com/AlexsLemonade/refinebio/issues/1487

from django.core.management.base import BaseCommand
from nomad import Nomad
from django.db.models import Count, OuterRef, Subquery
from django.db.models.expressions import Q

from data_refinery_common.models import (
    Experiment,
    ProcessorJob,
    Sample,
)
from data_refinery_common.job_lookup import ProcessorEnum, ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.rna_seq import get_quant_results_for_experiment
from data_refinery_common.job_management import create_processor_job_for_original_files
from data_refinery_foreman.foreman import main

logger = get_and_configure_logger(__name__)

def requeue_downloader_jobs():
    """Creates a tximport job for all eligible experiments."""
    latest_processor_job_for_sample = ProcessorJob.objects\
            .filter(origina_files__samples=OuterRef('id'))\
            .order_by('-start_time')

    eligible_sample_ids = Sample.objects\
      .annotate(
        failure_reason=Subquery(latest_processor_job_for_sample.values('failure_reason')[:1])
      )\
      .filter(
        failure_reason__startswith='ProcessorJob has already completed'
      )\
      .values_list('id', flat=True)

    eligible_samples = Sample.objects.filter(id__in=eligible_sample_ids)\
                             .prefetch_related('original_files')

    total_samples_queued = 0
    for sample in eligible_samples:
      original_files = list(sample.original_files.all())
      if not len(original_files): continue

      create_downloader_job(original_files)
      total_samples_queued += 1

    logger.info("Re-queued %d samples that had failed .", total_samples_queued)

class Command(BaseCommand):
  def handle(self, *args, **options):
    requeue_downloader_jobs()
