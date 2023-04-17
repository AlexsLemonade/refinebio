"""
This command will import metadata that has been supplied by an external
contributor in our JSON format and apply it to our samples. NOTE: It assumes
that the metadata is properly formatted, so for some sources we might have to
write some glue code to put it into the correct format before importing.

The JSON format is documented in
https://github.com/AlexsLemonade/refinebio/issues/2127#issuecomment-591651893

This management command strives to be idempotent, so it should be safe to run
it again on the same input *as long as the submitter name remains the same*.
This means we can also use it to update our samples if someone has updated
metadata for us.
"""
import json
import sys
import uuid
from typing import Dict, List

from django.core.management.base import BaseCommand

import boto3
import botocore

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import *
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)


def import_sample_attributes(accession_code: str, attributes: List, source: Contribution):
    try:
        sample = Sample.objects.get(accession_code=accession_code)
    except Sample.DoesNotExist:
        logger.debug(
            "Skipping metadata for sample {} that we don't know about".format(accession_code)
        )
        return

    for attribute in attributes:
        if type(attribute) != dict:
            logger.error(
                "An observation for sample '{}' is of type '{}', not dict. Skipping...".format(
                    accession_code, type(attribute)
                )
            )
            continue

        for name, value in attribute.items():
            # We want importing sample attributes to be idempotent, so if an
            # attribute with the same submitter and name already exists, let's
            # update that instead.
            #
            # NOTE: Are there any situations where a sample could have multiple
            # attributes with the same ontology term name?
            try:
                attribute = SampleAttribute.objects.get(
                    sample=sample, source=source, name__ontology_term=name
                )
            except SampleAttribute.DoesNotExist:
                attribute = SampleAttribute()
                attribute.sample = sample
                attribute.source = source
                attribute.name = OntologyTerm.get_or_create_from_api(name)

            attribute.set_value(value["value"])

            probability = value.get("probability", None)
            if probability is not None:
                attribute.probability = probability

            unit = value.get("unit", None)
            if unit is not None:
                attribute.unit = OntologyTerm.get_or_create_from_api(unit)

            attribute.save()


def import_metadata(metadata: Dict, source: Contribution):
    for sample in metadata:
        if type(sample["sample_accession"]) != str:
            logger.error(
                "The provided sample accession {} is not a string".format(
                    sample["sample_accession"]
                )
            )
            continue

        if type(sample["attributes"]) != list:
            # Don't print sample_attributes itself because it could be massive
            logger.error(
                "The provided attributes for sample {} is not a list".format(
                    sample["sample_accession"]
                )
            )
            continue

        import_sample_attributes(sample["sample_accession"], sample["attributes"], source)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help=(
                "A JSON file containing sample metadata to import in the correct format.\n"
                + "s3:// URLs are also accepted."
            ),
        )
        parser.add_argument("--source-name", type=str, help=("The name of the source"))
        parser.add_argument("--methods-url", type=str, help=("A link to this metadata's methods"))

    def handle(self, *args, **options):
        okay = True
        if options["file"] is None:
            logger.error("You must specify a file to import metadata from")
            okay = False
        if options["source_name"] is None:
            logger.error("You must specify a source name")
            okay = False
        if options["methods_url"] is None:
            logger.error("You must specify a methods url")
            okay = False
        if not okay:
            sys.exit(1)

        if "s3://" in options["file"]:
            bucket, key = parse_s3_url(options["file"])
            s3 = boto3.resource("s3")
            try:
                filepath = "/tmp/metadata_" + str(uuid.uuid4()) + ".json"
                s3.Bucket(bucket).download_file(key, filepath)
            except botocore.exceptions.ClientError as e:
                if e.response["Error"]["Code"] == "404":
                    logger.error("The remote file does not exist.")
                    raise
                else:
                    raise
        else:
            filepath = options["file"]

        with open(filepath) as file:
            metadata = json.load(file)

        if type(metadata) != list:
            logger.error("The provided metadata file is not a list of metadata information")
            sys.exit(1)

        source, _ = Contribution.objects.get_or_create(
            source_name=options["source_name"], methods_url=options["methods_url"]
        )

        import_metadata(metadata, source)
