"""
This command will import an entire ontology into our database so that we do not
hammer the Ontology Lookup Service with requests for a bunch of individual
ontology terms.
"""

import sys

from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OntologyTerm

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--ontology", type=str, help=("The prefix of the ontology to import (e.g. UO)."),
        )

    def handle(self, *args, **options):
        if options.get("ontology", None) is None:
            logger.error("You must provide an ontology term with --ontology")
            sys.exit(1)

        OntologyTerm.import_entire_ontology(options["ontology"])
