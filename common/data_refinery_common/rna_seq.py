from itertools import groupby
from typing import List

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputationalResult, Experiment, OrganismIndex
from data_refinery_common.utils import get_env_variable

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
EARLY_TXIMPORT_MIN_PERCENT = 0.80


def should_run_tximport(experiment: Experiment, results, is_tximport_job: bool):
    """ Returns whether or not the experiment is eligible to have tximport
    run on it.

    results is a set of ComputationalResults for the samples that had salmon quant run on them.
    """
    num_quantified = len(results)
    if num_quantified == 0:
        return False

    salmon_versions = set()
    for result in results:
        if result.organism_index.salmon_version:
            salmon_versions.add(result.organism_index.salmon_version)

    if len(salmon_versions) > 1:
        # Tximport requires that all samples are processed with the same salmon version
        # https://github.com/AlexsLemonade/refinebio/issues/1496
        return False

    eligible_samples = experiment.samples.filter(source_database="SRA", technology="RNA-SEQ")

    num_eligible_samples = eligible_samples.count()
    if num_eligible_samples == 0:
        return False

    percent_complete = num_quantified / num_eligible_samples

    if percent_complete == 1.0:
        # If an experiment is fully quantified then we should run
        # tximport regardless of its size.
        return True

    if (
        is_tximport_job
        and num_eligible_samples >= EARLY_TXIMPORT_MIN_SIZE
        and percent_complete >= EARLY_TXIMPORT_MIN_PERCENT
    ):
        return True
    else:
        return False


def get_quant_results_for_experiment(experiment: Experiment, filter_old_versions=True):
    """Returns a set of salmon quant results from `experiment`."""
    # Subquery to calculate quant results
    # https://docs.djangoproject.com/en/2.2/ref/models/expressions/#subquery-expressions
    all_results = ComputationalResult.objects.filter(sample__in=experiment.samples.all())

    if filter_old_versions:
        # Salmon version gets saved as what salmon outputs, which includes this prefix.
        current_salmon_version = "salmon " + get_env_variable("SALMON_VERSION", "0.13.1")
        organisms = experiment.organisms.all()
        organism_indices = OrganismIndex.objects.filter(
            salmon_version=current_salmon_version, organism__in=organisms
        )
        all_results.filter(organism_index__id__in=organism_indices.values("id"))

    all_results.select_related("computedfile").filter(
        computedfile__s3_bucket__isnull=False, computedfile__s3_key__isnull=False
    )

    def get_sample_id_set(result):
        return {sample.id for sample in result.samples.all()}

    latest_results = set()
    for k, group in groupby(sorted(list(all_results), get_sample_id_set), get_sample_id_set):
        latest_result = None
        for result in group:
            if not latest_result:
                latest_result = result
            else:
                if result.created_at > latest_result.created_at:
                    latest_result = result

        latest_results.add(latest_result)

    return latest_results


def get_quant_files_for_results(results: List[ComputationalResult]):
    """Returns a list of salmon quant results from `experiment`."""
    quant_files = []

    for result in results:
        quant_sf_file = result.get_quant_sf_file()
        if quant_sf_file:
            quant_files.append(quant_sf_file)
        else:
            logger.exception(
                "Salmon quant result found without quant.sf ComputedFile!",
                quant_result=result.id,
                sample=result.samples.first(),
            )
            raise Exception("Salmon quant result found without quant.sf ComputedFile!")

    return quant_files


ENA_DOWNLOAD_URL_TEMPLATE = (
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{short_accession}{sub_dir}"
    "/{long_accession}/{long_accession}{read_suffix}.fastq.gz"
)
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
        read_suffix=read_suffix,
    )
