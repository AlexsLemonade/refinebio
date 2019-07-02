from typing import Dict, List

from data_refinery_common.job_lookup import ProcessorEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Experiment,
)


def should_run_tximport(experiment: Experiment,
                        num_quantified: int,
                        is_tximport_job: bool):
    """ Returns whether or not the experiment is eligible to have tximport
    run on it.

    num_quantified is how many samples have had salmon quant run on them.
    """
    if num_quantified == 0:
        return False

    eligible_samples = experiment.samples.filter(source_database='SRA', technology='RNA-SEQ')

    num_eligible_samples = eligible_samples.count()
    if num_eligible_samples == 0:
        return False

    percent_complete = num_quantified / num_eligible_samples

    if percent_complete == 1.0:
        # If an experiment is fully quantified then we should run
        # tximport regardless of its size.
        return True

    if is_tximport_job \
       and num_eligible_samples >= EARLY_TXIMPORT_MIN_SIZE \
       and percent_complete >= EARLY_TXIMPORT_MIN_PERCENT:
        return True
    else:
        return False


def get_quant_results_for_experiment(experiment: Experiment):
    """Returns a list of salmon quant results from `experiment`."""
    results = []
    for sample in experiment.samples.all():
        for result in sample.results.order_by('-created_at').all():
            # TODO: this will break when we want to run for a new version.
            if result.processor.name == ProcessorEnum.SALMON_QUANT.value['name']:
                results.append(result)
                break

    return results


def get_quant_files_for_results(results: List[ComputationalResult]):
    """Returns a list of salmon quant results from `experiment`."""
    quant_files = []
    for result in results:
        try:
            quant_files.append(ComputedFile.objects.filter(
                result=result,
                filename="quant.sf",
                s3_key__isnull=False,
                s3_bucket__isnull=False,
                ).order_by('-id')[0])
        except Exception as e:
            try:
                sample = result.samples.first()
            except:
                sample = None

            logger.exception(
                "Salmon quant result found without quant.sf ComputedFile!",
                quant_result=result.id,
                sample=sample.id,
                experiment=experiment.id
            )
            raise e

    return quant_files
