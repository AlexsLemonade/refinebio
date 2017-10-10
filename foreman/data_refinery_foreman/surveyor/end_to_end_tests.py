import time
import json
from django.test import TransactionTestCase
from unittest.mock import patch, Mock
from data_refinery_common.models import (
    Organism,
    SurveyJob,
    SurveyJobKeyValue,
    Batch,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_foreman.surveyor import surveyor
from data_refinery_foreman.surveyor.array_express import EXPERIMENTS_URL

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


LOOP_TIME = 5  # seconds

EXPERIMENTS_JSON = """
{
  "experiments": {
    "experiment": [
      {
        "bioassaydatagroup": [
          {
            "isderived": 0,
            "bioassays": 1,
            "dataformat": "rawData",
            "arraydesignprovider": null,
            "bioassaydatacubes": 1,
            "name": "rawData",
            "id": null
          },
          {
            "isderived": 0,
            "bioassays": 1,
            "dataformat": "scan",
            "arraydesignprovider": null,
            "bioassaydatacubes": 1,
            "name": "scan",
            "id": null
          },
          {
            "isderived": 1,
            "bioassays": 1,
            "dataformat": "processedDataMatrix",
            "arraydesignprovider": null,
            "bioassaydatacubes": 1,
            "name": "processedDataMatrix",
            "id": null
          }
        ],
        "protocol": [
          {
            "accession": "P-GSE22166-6",
            "id": 52818
          },
          {
            "accession": "P-GSE22166-8",
            "id": 52820
          },
          {
            "accession": "P-GSE22166-5",
            "id": 52814
          },
          {
            "accession": "P-GSE22166-7",
            "id": 52819
          },
          {
            "accession": "P-GSE22166-4",
            "id": 52815
          },
          {
            "accession": "P-GSE22166-3",
            "id": 52816
          },
          {
            "accession": "P-GSE22166-2",
            "id": 52817
          },
          {
            "accession": "P-GSE22166-1",
            "id": 52813
          }
        ],
        "arraydesign": [
          {
            "legacy_id": 405156763,
            "count": 1,
            "name": "Affymetrix GeneChip Human Genome U133 Plus 2.0 [HG-U133_Plus_2]",
            "accession": "A-AFFY-44",
            "id": 11119
          }
        ],
        "samplecharacteristic": [
          {
            "value": [
              "Homo sapiens"
            ],
            "category": "Organism"
          },
          {
            "value": [
              "Human umbilical vein endothelial cell"
            ],
            "category": "tissue"
          }
        ],
        "provider": [
          {
            "email": null,
            "role": null,
            "contact": "Steven Krilis"
          },
          {
            "email": null,
            "role": null,
            "contact": "Pei Yu"
          },
          {
            "email": "j.ioannou@ich.ucl.ac.uk",
            "role": "submitter",
            "contact": "Yiannis Ioannou"
          }
        ],
        "description": [
          {
            "text": "The original text was too long. This is a placeholder.",
            "id": null
          }
        ],
        "experimenttype": [
          "transcription profiling by array"
        ],
        "organism": [
          "Homo sapiens"
        ],
        "lastupdatedate": "2011-06-10",
        "submissiondate": "2010-06-03",
        "releasedate": "2010-06-16",
        "name": "Oxidoreductase transcript data from human umbilical vein endothelial cells",
        "secondaryaccession": [
          "GSE22166"
        ],
        "accession": "E-GEOD-22166",
        "id": 16991
      }
    ],
    "total-assays": 1,
    "total-samples": 1,
    "total": 1,
    "revision": "091015",
    "version": 3.0,
    "api-revision": "091015",
    "api-version": 3
  }
} """

SAMPLES_JSON = """
{
  "experiment": {
    "sample": [
      {
        "source": {
          "comment": [
            {
              "value": "Gene expression data from embryos younger than.... truncated",
              "name": "Sample_description"
            },
            {
              "value": "Unstimulated HUVEC",
              "name": "Sample_source_name"
            }
          ],
          "name": "GSM551183 1"
        },
        "scan": {
          "name": "GSM551183.CEL"
        },
        "labeled-extract": {
          "label": "biotin",
          "name": "GSM551183 LE 1"
        },
        "file": [
          {
            "comment": {
              "value": "http://nginx/E-GEOD-22166.raw.1.zip",
              "name": "ArrayExpress FTP file"
            },
            "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-22166/E-GEOD-22166.raw.1.zip/GSM551183.CEL",
            "name": "GSM551183.CEL",
            "type": "data"
          },
          {
            "comment": {
              "value": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-22166/E-GEOD-22166.processed.1.zip",
              "name": "Derived ArrayExpress FTP file"
            },
            "url": "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-22166/E-GEOD-22166.processed.1.zip/GSM551183_sample_table.txt",
            "name": "GSM551183_sample_table.txt",
            "type": "derived data matrix"
          }
        ],
        "extract": {
          "name": "GSM551183 extract 1"
        },
        "characteristic": [
          {
            "value": "Homo sapiens",
            "category": "Organism"
          },
          {
            "value": "Human umbilical vein endothelial cell",
            "category": "tissue"
          }
        ],
        "assay": {
          "name": "GSM551183"
        }
      }
    ],
    "accession": "E-GEOD-22166",
    "revision": "091015",
    "version": 1.0,
    "api-revision": "091015",
    "api-version": 3
  }
} """


def mocked_requests_get(url):
    mock = Mock(ok=True)
    if url == (EXPERIMENTS_URL + "E-GEOD-22166"):
        mock.json.return_value = json.loads(EXPERIMENTS_JSON)
    else:
        mock.json.return_value = json.loads(SAMPLES_JSON)

    return mock


def wait_for_job(job, job_class: type):
    """Monitors the `job_class` table for when `job` is done."""
    job = job_class.objects.filter(id=job.id).get()
    while job.success is None:
        logger.info("Still polling the %s.",
                    job_class.__name__)
        time.sleep(LOOP_TIME)
        job = job_class.objects.filter(id=job.id).get()

    return job


# TransactionTestCase makes database calls complete before the test
# ends.  Otherwise the workers wouldn't actually be able to find the
# job in the database cause it'd be stuck in a transaction.
class ScanUpcEndToEndTestCase(TransactionTestCase):
    @patch('data_refinery_foreman.surveyor.array_express.requests.get')
    def test_calls_survey(self, mock_get):
        """If source_type is supported calls the appropriate survey method."""
        mock_get.side_effect = mocked_requests_get

        # Prevent a call being made to NCBI's API to determine
        # organism name/id.
        organism = Organism(name="HOMO SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value="E-GEOD-22166")
        key_value_pair.save()

        surveyor.run_job(survey_job)
        logger.info("Started Survey Job %d, waiting for it to complete.", survey_job.id)
        survey_job = wait_for_job(survey_job, SurveyJob)
        self.assertTrue(survey_job.success)

        batch = Batch.objects.all()[0]
        batch = Batch.objects.filter(survey_job=survey_job).get()

        downloader_job = batch.downloaderjob_set.get()
        logger.info("Survey Job finished, waiting for Downloader Job %d to complete.",
                    downloader_job.id)
        downloader_job = wait_for_job(downloader_job, DownloaderJob)
        self.assertTrue(downloader_job.success)

        processor_job = batch.processorjob_set.get()
        logger.info("Downloader Job finished, waiting for processor Job %d to complete.",
                    processor_job.id)
        processor_job = wait_for_job(processor_job, ProcessorJob)
        self.assertTrue(processor_job.success)
