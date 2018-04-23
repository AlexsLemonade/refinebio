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
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_foreman.surveyor.harmony import harmonize, parse_sdrf

class HarmonyTestCase(TestCase):
    def setUp(self):
        self.sample = Sample()
        self.sample.save()

        self.samples = [self.sample]

    def test_sdrf_harmony(self):
        """

        """
        
        metadata = parse_sdrf("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3050/E-MTAB-3050.sdrf.txt")
        for sample in metadata:
            sample['title'] = sample['Assay Name']
        harmonized = harmonize(metadata)
        print(harmonized)

    def test_sra_harmony(self):
        """

        """
        
        metadata = SraSurveyor.gather_all_metadata("SRR6718414")
        metadata['title'] = metadata['sample_title']
        harmonized = harmonize([metadata])
        print(harmonized)