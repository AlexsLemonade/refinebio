from django.test import TestCase
from data_refinery_common import file_management
from data_refinery_models.models import Batch, BatchStatuses, SurveyJob


class FilesTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        survey_job = SurveyJob(
            source_type="ARRAY_EXPRESS"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="ARRAY_EXPRESS",
            size_in_bytes=0,
            download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",  # noqa
            raw_format="CEL",
            processed_format="PCL",
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-141",
            experiment_accession_code="E-GEOD-59071",
            experiment_title="It doesn't really matter.",
            name="GSM1426072.CEL",
            internal_location="A-AFFY-141/AFFY_TO_PCL",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()

        super(FilesTestCase, cls).setUpClass()

    def test_raw_dir_s3(self):
        file_management.USE_S3 = True
        batch = Batch.objects.get(pk=1)
        raw_dir = file_management.get_raw_dir(batch)

        self.assertEqual(raw_dir, "raw/A-AFFY-141/AFFY_TO_PCL")

    def test_raw_dir_no_s3(self):
        file_management.USE_S3 = False
        batch = Batch.objects.get(pk=1)
        raw_dir = file_management.get_raw_dir(batch)

        self.assertEqual(raw_dir, "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL")

    def test_raw_download_path_s3(self):
        file_management.USE_S3 = True
        batch = Batch.objects.get(pk=1)
        download_path = file_management.get_raw_download_path(batch)
        correct_path = "raw/A-AFFY-141/AFFY_TO_PCL/E-GEOD-59071.raw.3.zip"

        self.assertEqual(download_path, correct_path)

    def test_raw_download_path_no_s3(self):
        file_management.USE_S3 = False
        batch = Batch.objects.get(pk=1)
        download_path = file_management.get_raw_download_path(batch)
        correct_path = "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/E-GEOD-59071.raw.3.zip"

        self.assertEqual(download_path, correct_path)

    def test_raw_path_s3(self):
        file_management.USE_S3 = True
        batch = Batch.objects.get(pk=1)
        raw_path = file_management.get_raw_path(batch)
        correct_path = "raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072.CEL"

        self.assertEqual(raw_path, correct_path)

    def test_raw_path_no_s3(self):
        file_management.USE_S3 = False
        batch = Batch.objects.get(pk=1)
        raw_path = file_management.get_raw_path(batch)
        correct_path = "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072.CEL"

        self.assertEqual(raw_path, correct_path)

    def test_temp_dir(self):
        batch = Batch.objects.get(pk=1)
        temp_dir = file_management.get_temp_dir(batch)

        self.assertEqual(temp_dir, "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1")

    def test_temp_pre_path(self):
        batch = Batch.objects.get(pk=1)
        temp_pre_path = file_management.get_temp_pre_path(batch)
        correct_path = "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1/GSM1426072.CEL"

        self.assertEqual(temp_pre_path, correct_path)

    def test_temp_post_path(self):
        batch = Batch.objects.get(pk=1)
        temp_post_path = file_management.get_temp_post_path(batch)
        correct_path = "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1/GSM1426072.PCL"

        self.assertEqual(temp_post_path, correct_path)

    def test_processed_dir_s3(self):
        file_management.USE_S3 = True
        batch = Batch.objects.get(pk=1)
        processed_dir = file_management.get_processed_dir(batch)

        self.assertEqual(processed_dir, "processed/A-AFFY-141/AFFY_TO_PCL")

    def test_processed_dir_no_s3(self):
        file_management.USE_S3 = False
        batch = Batch.objects.get(pk=1)
        processed_dir = file_management.get_processed_dir(batch)

        self.assertEqual(processed_dir, "/home/user/data_store/processed/A-AFFY-141/AFFY_TO_PCL")

    def test_processed_path_s3(self):
        file_management.USE_S3 = True
        batch = Batch.objects.get(pk=1)
        processed_path = file_management.get_processed_path(batch)
        correct_path = "processed/A-AFFY-141/AFFY_TO_PCL/GSM1426072.PCL"

        self.assertEqual(processed_path, correct_path)

    def test_processed_path_no_s3(self):
        file_management.USE_S3 = False
        batch = Batch.objects.get(pk=1)
        processed_path = file_management.get_processed_path(batch)
        correct_path = "/home/user/data_store/processed/A-AFFY-141/AFFY_TO_PCL/GSM1426072.PCL"

        self.assertEqual(processed_path, correct_path)
