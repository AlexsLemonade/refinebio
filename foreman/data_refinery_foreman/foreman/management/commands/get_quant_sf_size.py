import csv

from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Experiment

logger = get_and_configure_logger(__name__)

csv_headers = ["sample_accession_code", "experiment_accession_code", "organism", "size_in_bytes"]


class Command(BaseCommand):
    def handle(self, *args, **options):
        with open("/tmp/quant_sf_file_sizes.csv", "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_headers)
            writer.writeheader()

            rnaseq_experiments = Experiment.objects.filter(technology="RNA-SEQ").prefetch_related(
                "samples__results"
            )
            for experiment in rnaseq_experiments:
                logger.info("Crunching experiment %s", experiment.accession_code)
                for sample in experiment.samples.all():
                    # Some experiments are multi-technology, but we only want RNASeq samples.
                    if sample.technology == "RNA-SEQ":
                        for result in sample.results.all():
                            quant_sf_file = result.get_quant_sf_file()
                            if quant_sf_file:
                                writer.writerow(
                                    {
                                        "sample_accession_code": sample.accession_code,
                                        "experiment_accession_code": experiment.accession_code,
                                        "organism": sample.organism.name,
                                        "size_in_bytes": quant_sf_file.size_in_bytes,
                                    }
                                )
