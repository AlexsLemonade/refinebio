import json
import datetime
from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism,
    Sample
)
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor
from data_refinery_foreman.surveyor.harmony import harmonize, parse_sdrf

class HarmonyTestCase(TestCase):
    def setUp(self):
        self.sample = Sample()
        self.sample.save()

        self.samples = [self.sample]

    def test_harmony(self):
        """

        """
        
        metadata = parse_sdrf("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3050/E-MTAB-3050.sdrf.txt")
        harmonized = harmonize(metadata)
        print(harmonized)