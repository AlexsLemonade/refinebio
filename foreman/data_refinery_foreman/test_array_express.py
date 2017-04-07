import json
from unittest.mock import Mock, patch
from django.test import TestCase
from data_refinery_models.models import (
    Batch,
    SurveyJob,
    SurveyJobKeyValue
)
import data_refinery_foreman.surveyor as surveyor
from data_refinery_foreman.surveyor.array_express_surveyor \
    import ArrayExpressSurveyor


class SurveyTestCase(TestCase):
    experiments_json = """
    {
    "files": {
        "api-revision": "091015",
        "api-version": 2,
        "experiment": [
            {
                "accession": "E-MTAB-3050",
                "file": [
                    {
                        "extension": "zip",
                        "kind": "raw",
                        "lastmodified": "2014-10-30T10:15:00",
                        "location": "E-MTAB-3050.raw.1.zip",
                        "name": "E-MTAB-3050.raw.1.zip",
                        "size": 14876114,
                        "url":
    "http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3050/E-MTAB-3050.raw.1.zip"
                    },
                    {
                        "extension": "xls",
                        "kind": "adf",
                        "lastmodified": "2010-03-14T02:31:00",
                        "location": "A-AFFY-1.adf.xls",
                        "name": "A-AFFY-1.adf.xls",
                        "size": 2040084,
                        "url":
    "http://www.ebi.ac.uk/arrayexpress/files/A-AFFY-1/A-AFFY-1.adf.xls"
                    }
                ]
            },
            {
                "accession": "E-MTAB-3042",
                "file": [
                    {
                        "extension": "txt",
                        "kind": "idf",
                        "lastmodified": "2014-10-28T10:15:00",
                        "location": "E-MTAB-3042.idf.txt",
                        "name": "E-MTAB-3042.idf.txt",
                        "size": 5874,
                        "url":
    "http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3042/E-MTAB-3042.idf.txt"
                    },
                    {
                        "extension": "zip",
                        "kind": "raw",
                        "lastmodified": "2014-10-28T10:15:00",
                        "location": "E-MTAB-3042.raw.1.zip",
                        "name": "E-MTAB-3042.raw.1.zip",
                        "size": 5525709,
                        "url":
    "http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3042/E-MTAB-3042.raw.1.zip"
                    }
                ]
            }
        ],
        "revision": 130311,
        "total-experiments": 108,
        "version": 1.2
    }
}
"""

    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="accession_code",
                                           value="A-AFFY-1")
        key_value_pair.save()

    def tearDown(self):
        SurveyJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        Batch.objects.all().delete

    @patch('data_refinery_foreman.surveyor.array_express_surveyor.requests.get')
    def test_multiple_experiements(self, mock_get):
        """Multiple experiments are turned into multiple batches"""
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.json.return_value = json.loads(
            self.experiments_json)

        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        self.assertTrue(ae_surveyor.survey(self.survey_job))
        self.assertEqual(2, Batch.objects.all().count())

    experiment_json = """
    {
    "files": {
        "api-revision": "091015",
        "api-version": 2,
        "experiment": [
            {
                "accession": "E-MTAB-3050",
                "file": {
                        "extension": "zip",
                        "kind": "raw",
                        "lastmodified": "2014-10-30T10:15:00",
                        "location": "E-MTAB-3050.raw.1.zip",
                        "name": "E-MTAB-3050.raw.1.zip",
                        "size": 14876114,
                        "url":
        "http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3050/E-MTAB-3050.raw.1.zip"
                    }
                }
    ],
        "revision": 130311,
        "total-experiments": 108,
        "version": 1.2
    }
}
"""

    @patch('data_refinery_foreman.surveyor.array_express_surveyor.requests.get')
    def test_multiple_experiements(self, mock_get):
        """A single experiment is turned into a single batch."""
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.json.return_value = json.loads(
            self.experiment_json)

        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        self.assertTrue(ae_surveyor.survey(self.survey_job))
        self.assertEqual(1, Batch.objects.all().count())
