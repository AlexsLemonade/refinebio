from typing import Dict, List

from data_refinery_common.job_lookup import ProcessorEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Experiment,
)


logger = get_and_configure_logger(__name__)


# Some experiments won't be entirely processed, but we'd still like to
# make the samples we can process available. This means we need to run
# tximport on the experiment before 100% of the samples are processed
# individually.
# This idea has been discussed here: https://github.com/AlexsLemonade/refinebio/issues/909

# The consensus is that this is a good idea, but that we need a cutoff
# to determine which experiments have enough data to have tximport run
# on them early.  Candace ran an experiment to find these cutoff
# values and recorded the results of this experiment here:
# https://github.com/AlexsLemonade/tximport_partial_run_tests/pull/3

# The gist of that discussion/experiment is that we need two cutoff
# values, one for a minimum size experiment that can be processed
# early and the percentage of completion necessary before we can
# run tximport on the experiment. The values we decided on are:
EARLY_TXIMPORT_MIN_SIZE = 25
EARLY_TXIMPORT_MIN_PERCENT = .80


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


ENA_DOWNLOAD_URL_TEMPLATE = ("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{short_accession}{sub_dir}"
                             "/{long_accession}/{long_accession}{read_suffix}.fastq.gz")
ENA_SUB_DIR_PREFIX = "/00"


def _build_ena_file_url(run_accession: str, read_suffix=""):
    # ENA has a weird way of nesting data: if the run accession is
    # greater than 9 characters long then there is an extra
    # sub-directory in the path which is "00" + the last digit of
    # the run accession.
    sub_dir = ""
    if len(run_accession) > 9:
        sub_dir = ENA_SUB_DIR_PREFIX + run_accession[-1]

    return ENA_DOWNLOAD_URL_TEMPLATE.format(
        short_accession=run_accession[:6],
        sub_dir=sub_dir,
        long_accession=run_accession,
        read_suffix=read_suffix)
