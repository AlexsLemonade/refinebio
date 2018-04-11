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
from data_refinery_foreman.surveyor.array_express import (
    ArrayExpressSurveyor,
    EXPERIMENTS_URL,
)


EXPERIMENTS_JSON = """
{
    "experiments": {
        "api-revision": "091015",
        "api-version": 3,
        "experiment": [
            {
                "accession": "E-MTAB-3050",
                "arraydesign": [
                    {
                        "accession": "A-AFFY-1",
                        "count": 5,
                        "id": 11048,
                        "legacy_id": 5728564,
                        "name": "Affymetrix GeneChip Human  U95Av2 [HG_U95Av2]"
                    }
                ],
                "bioassaydatagroup": [
                    {
                        "arraydesignprovider": null,
                        "bioassaydatacubes": 5,
                        "bioassays": 5,
                        "dataformat": "rawData",
                        "id": null,
                        "isderived": 0,
                        "name": "rawData"
                    }
                ],
                "description": [
                    {
                        "id": null,
                        "text": "description tex"
                    }
                ],
                "experimentalvariable": [
                    {
                        "name": "cell type",
                        "value": [
                            "differentiated",
                            "expanded",
                            "freshly isolated"
                        ]
                    }
                ],
                "experimentdesign": [
                    "cell type comparison design",
                    "development or differentiation design"
                ],
                "experimenttype": [
                    "transcription profiling by array"
                ],
                "id": 511696,
                "lastupdatedate": "2014-10-30",
                "name": "Microarray analysis of in vitro differentiation",
                "organism": [
                    "Homo sapiens"
                ],
                "protocol": [
                    {
                        "accession": "P-MTAB-41859",
                        "id": 1092859
                    }
                ],
                "provider": [
                    {
                        "contact": "Joel Habener",
                        "email": "jhabener@partners.org",
                        "role": "submitter"
                    }
                ],
                "releasedate": "2014-10-31",
                "samplecharacteristic": [
                    {
                        "category": "age",
                        "value": [
                            "38 year",
                            "54 year"
                        ]
                    }
                ]
            }
        ],
        "revision": "091015",
        "total": 1,
        "total-assays": 5,
        "total-samples": 2,
        "version": 3.0
    }
} """

SAMPLES_JSON = """
{
    "experiment": {
        "accession": "E-MTAB-3050",
        "api-revision": "091015",
        "api-version": 3,
        "revision": "091015",
        "sample": [
            {
                "assay": {
                    "name": "1007409-C30057"
                },
                "characteristic": [
                    {
                        "category": "organism",
                        "value": "Homo sapiens"
                    }
                ],
                "extract": {
                    "name": "donor A islets RNA"
                },
                "file": [
                    {
                        "comment": {
                            "name": "ArrayExpress FTP file",
                            "value": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip"
                        },
                        "name": "C30057.CEL",
                        "type": "data",
                        "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip/C30057.CEL"
                    },
                    {
                        "comment": {
                            "name": "Derived ArrayExpress FTP file",
                            "value": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.processed.1.zip"
                        },
                        "name": "C30057.txt",
                        "type": "derived data",
                        "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.processed.1.zip/C30057.txt"
                    }
                ],
                "labeled-extract": {
                    "label": "biotin",
                    "name": "donor A islets LEX"
                },
                "source": {
                    "name": "donor A islets"
                },
                "variable": [
                    {
                        "name": "cell type",
                        "value": "freshly isolated"
                    }
                ]
            },
            {
                "assay": {
                    "name": "1007409-C30058"
                },
                "characteristic": [
                    {
                        "category": "organism",
                        "value": "Homo sapiens"
                    }
                ],
                "extract": {
                    "name": "donor A expanded cells RNA"
                },
                "file": [
                    {
                        "comment": {
                            "name": "ArrayExpress FTP file",
                            "value": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip"
                        },
                        "name": "C30058.CEL",
                        "type": "data",
                        "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip/C30058.CEL"
                    },
                    {
                        "comment": {
                            "name": "Derived ArrayExpress FTP file",
                            "value": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.processed.1.zip"
                        },
                        "name": "C30058.txt",
                        "type": "derived data",
                        "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.processed.1.zip/C30058.txt"
                    }
                ],
                "labeled-extract": {
                    "label": "biotin",
                    "name": "donor A expanded cells LEX"
                },
                "source": {
                    "name": "donor A islets"
                },
                "variable": [
                    {
                        "name": "cell type",
                        "value": "expanded"
                    }
                ]
            }
        ],
        "version": 1.0
    }
}"""  # noqa


def mocked_requests_get(url):
    mock = Mock(ok=True)
    if url == (EXPERIMENTS_URL + "E-MTAB-3050"):
        mock.json.return_value = json.loads(EXPERIMENTS_JSON)
    else:
        mock.json.return_value = json.loads(SAMPLES_JSON)

    return mock


class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value="E-MTAB-3050")
        key_value_pair.save()

        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(name="HOMO_SAPIENS",
                            taxonomy_id=9606,
                            is_scientific_name=True)
        organism.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    # @patch('data_refinery_foreman.surveyor.array_express.requests.get')
    @patch('data_refinery_foreman.surveyor.external_source.send_job')
    def test_survey(self, mock_send_task):
        """The 'survey' function generates one Batch per sample.

        This test also tests the handle_batches method of ExternalSourceSurveyor
        which isn't tested on its own because it is an abstract class.
        """
        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        ae_surveyor.survey()

        samples = Sample.objects.all()
        downloader_jobs = DownloaderJob.objects.all()

        # We are expecting this to discoever 5 samples.
        self.assertEqual(samples.count(), 5)
        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 1)
