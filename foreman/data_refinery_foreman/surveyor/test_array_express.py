import json
import datetime
from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    Batch,
    File,
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism
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
        organism = Organism(name="HOMO SAPIENS",
                            taxonomy_id=9606,
                            is_scientific_name=True)
        organism.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        Batch.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    @patch('data_refinery_foreman.surveyor.array_express.requests.get')
    def test_experiment_object(self, mock_get):
        """The get_experiment_metadata function extracts all experiment metadata
        from the experiments API."""
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.json.return_value = json.loads(EXPERIMENTS_JSON)

        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        experiment = ae_surveyor.get_experiment_metadata("E-MTAB-3050")
        self.assertEqual("Microarray analysis of in vitro differentiation", experiment["name"])
        self.assertEqual("E-MTAB-3050", experiment["experiment_accession_code"])
        self.assertEqual("A-AFFY-1", experiment["platform_accession_code"])
        self.assertEqual("2014-10-31", experiment["release_date"])
        self.assertEqual("2014-10-30", experiment["last_update_date"])

    @patch('data_refinery_foreman.surveyor.array_express.requests.get')
    @patch('data_refinery_foreman.surveyor.external_source.send_job')
    def test_survey(self, mock_send_task, mock_get):
        """The 'survey' function generates one Batch per sample.

        This test also tests the handle_batches method of ExternalSourceSurveyor
        which isn't tested on its own because it is an abstract class.
        """
        mock_send_task.return_value = Mock(ok=True)
        mock_get.side_effect = mocked_requests_get

        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        ae_surveyor.survey()

        downloader_jobs = DownloaderJob.objects.all()
        mock_send_task.assert_has_calls([
            call(Downloaders.ARRAY_EXPRESS,
                 downloader_jobs[0].id)
        ])
        batches = Batch.objects.all()
        self.assertEqual(2, len(batches))
        self.assertEqual(1, len(downloader_jobs))

        batch = batches[0]
        self.assertEqual(batch.survey_job.id, self.survey_job.id)
        self.assertEqual(batch.source_type, "ARRAY_EXPRESS")
        self.assertEqual(batch.pipeline_required, "AFFY_TO_PCL")
        self.assertEqual(batch.platform_accession_code, "A-AFFY-1")
        self.assertEqual(batch.experiment_accession_code, "E-MTAB-3050")
        self.assertEqual(batch.experiment_title, "Microarray analysis of in vitro differentiation")
        self.assertEqual(batch.status, "NEW")
        self.assertEqual(batch.release_date, datetime.date(2014, 10, 31))
        self.assertEqual(batch.last_uploaded_date, datetime.date(2014, 10, 30))
        self.assertEqual(batch.organism_id, 9606)
        self.assertEqual(batch.organism_name, "HOMO SAPIENS")

        file = batch.files[0]
        self.assertEqual(file.size_in_bytes, -1)
        self.assertEqual(file.download_url, "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip")  # noqa
        self.assertEqual(file.raw_format, "CEL")
        self.assertEqual(file.processed_format, "PCL")
        self.assertEqual(file.name, "C30057.CEL")
        self.assertEqual(file.internal_location, "A-AFFY-1/AFFY_TO_PCL")


class GroupBatchesTestCase(TestCase):
    def test_survey(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        surveyor = ArrayExpressSurveyor(survey_job)
        file1 = File(download_url="a")
        file2 = File(download_url="a")
        file3 = File(download_url="b")
        file4 = File(download_url="a")
        batch1 = Batch(files=[file1])
        batch2 = Batch(files=[file2])
        batch3 = Batch(files=[file3, file4])

        surveyor.batches = [batch1, batch2, batch3]
        groups = surveyor.group_batches()
        self.assertEqual(groups, [[batch1, batch2], [batch3]])
