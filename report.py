import os
import numpy as np
import pandas as pd
import shutil
from django.db.models import Count
from data_refinery_common.models import Organism, Sample
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)
ROOT_DIR = "/home/user/data_store/"


def organism_can_have_qn_target(organism: Organism, sample_threshold=100):
    """ Check that the organism has more than `sample_threshold` samples on
    some microarray platform """
    microarray_platforms = organism.sample_set\
        .filter(has_raw=True, technology="MICROARRAY", is_processed=True)\
        .values('platform_accession_code')\
        .annotate(count=Count('id'))\
        .filter(count__gt=sample_threshold)

    return microarray_platforms.exists()


organism_name = "DANIO_RERIO"
if organism_name:
    organisms = [Organism.get_object_for_name(organism_name.upper())]
else:
    organisms = Organism.objects.exclude(name__in=['HOMO_SAPIENS', 'MUS_MUSCULUS'])

for organism in organisms:
    if not organism_can_have_qn_target(organism):
        continue

    organism_dir = os.path.join(ROOT_DIR, organism.name)
    os.makedirs(organism_dir, exist_ok=True)

    samples = organism.sample_set.filter(has_raw=True, technology="MICROARRAY", is_processed=True)
    if samples.count() == 0:
        logger.error("No processed samples for organism.",
                     organism=organism,
                     count=samples.count())
        continue

    platform_counts = samples.values(
        'platform_accession_code'
    ).annotate(
        dcount=Count('platform_accession_code')
    ).order_by('-dcount')
    biggest_platform = platform_counts[0]['platform_accession_code']

    platform_samples = Sample.processed_objects.filter(
        platform_accession_code=biggest_platform,
        has_raw=True,
        technology="MICROARRAY",
        organism=organism,
        is_processed=True
    )

    report_filename = ROOT_DIR + organism.name + "_NAs.csv"
    with open(report_filename, "w") as report_file:
        report_file.write("sample_accession_code,num_NA_values\n")
        for sample in platform_samples:
            try:
                computed_file = sample.get_most_recent_smashable_result_file(sample)
                path = computed_file.sync_from_s3(os.path.join(organism.name,
                                                               computed_file.filename))
                data = pd.read_csv(path,
                                   sep='\t',
                                   header=0,
                                   index_col=0,
                                   dtype={0: str, 1: np.float32},
                                   error_bad_lines=False)

                num_nans = data.isnull().sum(axis=1)[1]

                report_file.write("{},{}\n".format(sample.accession_code, num_nans))
                os.remove(path)

            except Exception as e:
                logger.exception("Sample %s caused an exception, continuing.",
                                 sample.accession_code)
                continue

        shutil.rmtree(organism_dir, ignore_errors=True)
