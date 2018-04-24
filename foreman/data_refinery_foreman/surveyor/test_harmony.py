import datetime
import GEOparse
import json

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
from data_refinery_foreman.surveyor.geo import GeoSurveyor
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

        self.assertTrue('1007409-C30057' in harmonized.keys())
        self.assertTrue('sex' in harmonized['1007409-C30057'].keys())
        self.assertTrue('female' == harmonized['1007409-C30057']['sex'])
        self.assertTrue('age' in harmonized['1007409-C30057'].keys())
        self.assertTrue(54 == harmonized['1007409-C30057']['age'])
        self.assertTrue('part' in harmonized['1007409-C30057'].keys())
        self.assertTrue('subject' in harmonized['1007409-C30057'].keys())
        self.assertTrue('developmental_stage' in harmonized['1007409-C30057'].keys())

    def test_sra_harmony(self):
        """

        """
        
        metadata = SraSurveyor.gather_all_metadata("SRR6718414")
        metadata['title'] = metadata['sample_title']
        harmonized = harmonize([metadata])
        print(harmonized)

    def test_geo_harmony(self):
        """

        """
        
        # Illumina
        gse = GEOparse.get_GEO("GSE32628", destdir='/tmp')

        preprocessed_samples = []
        for sample_id, sample in gse.gsms.items():
            new_sample = {}
            for key, value in sample.metadata.items():

                if key == "characteristics_ch1":
                    for pair in value:
                        split = pair.split(':')
                        new_sample[split[0].strip()] = split[1].strip()
                    continue

                new_sample[key] = value[0]
            preprocessed_samples.append(new_sample)
        harmonized = harmonize(preprocessed_samples)
        

        self.assertTrue('SCC_P-57' in harmonized.keys())
        self.assertTrue('sex' in harmonized['SCC_P-57'].keys())
        self.assertTrue('age' in harmonized['SCC_P-57'].keys())
        self.assertTrue('part' in harmonized['SCC_P-57'].keys())
        self.assertTrue('subject' in harmonized['SCC_P-57'].keys())

        # Agilent Two Color 
        gse = GEOparse.get_GEO("GSE93857", destdir='/tmp')
        preprocessed_samples = []
        for sample_id, sample in gse.gsms.items():
            new_sample = {}
            for key, value in sample.metadata.items():

                if key == "characteristics_ch1":
                    for pair in value:
                        split = pair.split(':')
                        new_sample[split[0].strip()] = split[1].strip()
                    continue

                new_sample[key] = value[0]
            preprocessed_samples.append(new_sample)
        harmonized = harmonize(preprocessed_samples)
    

        gse = GEOparse.get_GEO("GSE103060", destdir='/tmp')
        preprocessed_samples = []
        for sample_id, sample in gse.gsms.items():
            new_sample = {}
            for key, value in sample.metadata.items():

                if key == "characteristics_ch1":
                    for pair in value:
                        split = pair.split(':')
                        new_sample[split[0].strip()] = split[1].strip()
                    continue

                new_sample[key] = value[0]
            preprocessed_samples.append(new_sample)
        harmonized = harmonize(preprocessed_samples)
        print(harmonized)
