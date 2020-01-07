from typing import Dict, List

from django.db.models import OuterRef, Subquery

from data_refinery_common.job_lookup import ProcessorEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputationalResult, ComputedFile, Experiment
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

    results is a queryset of ComputationalResults for the samples that had salmon quant run on them.
    """
    num_quantified = results.count()
    if num_quantified == 0:
        return False

    num_salmon_versions = (
        results.filter(organism_index__salmon_version__isnull=False)
        .values_list("organism_index__salmon_version")
        .distinct()
        .count()
    )
    if num_salmon_versions > 1:
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
    """Returns a list of salmon quant results from `experiment`."""
    # Subquery to calculate quant results
    # https://docs.djangoproject.com/en/2.2/ref/models/expressions/#subquery-expressions

    # Salmon version gets saved as what salmon outputs, which includes this prefix.
    current_salmon_version = "salmon " + get_env_variable("SALMON_VERSION", "0.13.1")

    if filter_old_versions:
        eligible_results = ComputationalResult.objects.prefetch_related("organism_index").filter(
            organism_index__salmon_version=current_salmon_version
        )
    else:
        eligible_results = ComputationalResult.objects.all()

    # A result is only eligible to be used if it actually got uploaded.
    eligible_results = eligible_results.select_related("computedfile").filter(
        computedfile__s3_bucket__isnull=False, computedfile__s3_key__isnull=False
    )

    # Calculate the computational results sorted that are associated with a given sample (
    # referenced from the top query)
    newest_computational_results = eligible_results.filter(
        samples=OuterRef("id"), processor__name=ProcessorEnum.SALMON_QUANT.value["name"],
    ).order_by("-created_at")

    # Annotate each sample in the experiment with the id of the most recent computational result
    computational_results_ids = (
        experiment.samples.all()
        .annotate(
            latest_computational_result_id=Subquery(newest_computational_results.values("id")[:1])
        )
        .filter(latest_computational_result_id__isnull=False)
        .values_list("latest_computational_result_id", flat=True)
    )

    # return the computational results that match those ids
    return ComputationalResult.objects.all().filter(id__in=computational_results_ids)


def get_quant_files_for_results(results: List[ComputationalResult]):
    """Returns a list of salmon quant results from `experiment`."""
    quant_files = []

    for result in results:
        try:
            quant_files.append(
                ComputedFile.objects.filter(
                    result=result,
                    filename="quant.sf",
                    s3_key__isnull=False,
                    s3_bucket__isnull=False,
                ).order_by("-id")[0]
            )
        except Exception as e:
            try:
                sample = result.samples.first()
            except:
                sample = None

            logger.exception(
                "Salmon quant result found without quant.sf ComputedFile!",
                quant_result=result.id,
                sample=sample.id,
                experiment=experiment.id,
            )
            raise e

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
