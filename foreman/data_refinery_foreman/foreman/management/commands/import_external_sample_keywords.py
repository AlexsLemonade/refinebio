"""
This command will import keywords that have been supplied by an external
contributor in our JSON format and apply it to our samples. NOTE: It assumes
that the input file is properly formatted, so for some sources we might have to
write some glue code to put it into the correct format before importing.

The JSON format is as follows:
    {
        "<ACCESSION>": ["<KEYWORD>", ...],
        ...
    }

This management command is idempotent, so it should be safe to run it again on
the same input *as long as the submitter name remains the same*. This means we
can also use it to update our samples if someone has updated keywords for us.
"""
import json
import sys
import uuid
from typing import Dict

from django.core.management.base import BaseCommand

import boto3
import botocore

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import *
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)


def import_keywords(keywords: Dict, source: Contribution):
    dirty_experiments = set()

    for accession_code, keywords in keywords.items():
        if type(accession_code) != str:
            logger.error("The provided sample accession {} is not a string".format(accession_code))
            continue

        if type(keywords) != list:
            # Don't print sample_attributes itself because it could be massive
            logger.error("The provided keywords for sample {} is not a list".format(accession_code))
            continue

        try:
            sample = Sample.objects.get(accession_code=accession_code)
        except Sample.DoesNotExist:
            logger.debug(
                "Skipping metadata for sample {} that we don't know about".format(accession_code)
            )
            continue

        dirty_experiments |= set(sample.experiment_set.all())

        for keyword in keywords:
            SampleKeyword.objects.get_or_create(
                sample=sample,
                source=source,
                name=OntologyTerm.get_or_create_from_api(keyword),
            )

    for experiment in dirty_experiments:
        experiment.update_sample_metadata_fields()
        experiment.save()


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help=(
                "A JSON file containing sample keywords to import in the correct format.\n"
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
                filepath = "/tmp/keyword_" + str(uuid.uuid4()) + ".json"
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
            keywords = json.load(file)

        if type(keywords) != dict:
            logger.error("The provided keywords file is not a dict with accession code keys")
            sys.exit(1)

        source, _ = Contribution.objects.get_or_create(
            source_name=options["source_name"], methods_url=options["methods_url"]
        )

        import_keywords(keywords, source)
