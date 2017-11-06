from __future__ import absolute_import, unicode_literals
import os
import string
import warnings
from typing import Dict
import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError
from celery import shared_task
from data_refinery_common.models import File
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


GENE_TO_TRANSCRIPT_TEMPLATE = "{gene_id}\t{transcript_id}\n"
GENE_TYPE_COLUMN = 2
# Removes each occurrance of ; and "
IDS_CLEANUP_TABLE = str.maketrans({";": None, "\"": None})


def _process_gtf(input_gtf_path: str,
                 output_gtf_path: str,
                 genes_to_transcripts_path: str) -> None:
    """Reads in a .gtf file and generates two new files from it.

    The first is a new .gtf file which has all of the pseudogenes
    filtered out of it. The other is a tsv mapping between gene_ids
    and transcript_ids.
    """
    with open(input_gtf_path, 'r') as input_gtf:
        with open(output_gtf_path, "w") as output_gtf:
            with open(genes_to_transcripts_path, "w") as genes_to_transcripts:
                for line in input_gtf:
                    # Filter out any lines containing "pseudogene".
                    if "pseudogene" in line:
                        continue

                    output_gtf.write(line)
                    tab_split_line = line.split("\t")

                    # Skip header lines:
                    if len(tab_split_line) < 2:
                        continue

                    if tab_split_line[GENE_TYPE_COLUMN] == "transcript":
                        ids_column = tab_split_line[-1].translate(IDS_CLEANUP_TABLE)
                        split_ids_column = ids_column.split(" ")

                        gene_id_index = split_ids_column.index("gene_id") + 1
                        gene_id = split_ids_column[gene_id_index]
                        transcript_id_index = split_ids_column.index("transcript_id") + 1
                        transcript_id = split_ids_column[transcript_id_index]

                        genes_to_transcripts.write(
                            GENE_TO_TRANSCRIPT_TEMPLATE.format(gene_id=gene_id,
                                                               transcript_id=transcript_id))
