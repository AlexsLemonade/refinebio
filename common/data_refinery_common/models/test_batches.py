from django.test import TestCase
from data_refinery_common.models import batches
from data_refinery_common.models import Batch, BatchStatuses, File, SurveyJob


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
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-141",
            experiment_accession_code="E-GEOD-59071",
            experiment_title="It doesn't really matter.",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()

        file = File(
            size_in_bytes=0,
            download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",  # noqa
            raw_format="CEL",
            processed_format="PCL",
            name="GSM1426072.CEL",
            internal_location="A-AFFY-141/AFFY_TO_PCL",
            batch=batch
        )
        file.save()

        super(FilesTestCase, cls).setUpClass()

    def test_raw_dir_s3(self):
        batches.USE_S3 = True
        file = File.objects.get(batch__id=1)
        raw_dir = file.get_raw_dir()

        self.assertEqual(raw_dir, "raw/A-AFFY-141/AFFY_TO_PCL")

    def test_raw_dir_no_s3(self):
        batches.USE_S3 = False
        file = File.objects.get(batch__id=1)
        raw_dir = file.get_raw_dir()

        self.assertEqual(raw_dir, "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL")

    def test_raw_download_path_s3(self):
        batches.USE_S3 = True
        file = File.objects.get(batch__id=1)
        download_path = file.get_raw_download_path()
        correct_path = "raw/A-AFFY-141/AFFY_TO_PCL/E-GEOD-59071.raw.3.zip"

        self.assertEqual(download_path, correct_path)

    def test_raw_download_path_no_s3(self):
        batches.USE_S3 = False
        file = File.objects.get(batch__id=1)
        download_path = file.get_raw_download_path()
        correct_path = "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/E-GEOD-59071.raw.3.zip"

        self.assertEqual(download_path, correct_path)

    def test_raw_path_s3(self):
        batches.USE_S3 = True
        file = File.objects.get(batch__id=1)
        raw_path = file.get_raw_path()
        correct_path = "raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072.CEL"

        self.assertEqual(raw_path, correct_path)

    def test_raw_path_no_s3(self):
        batches.USE_S3 = False
        file = File.objects.get(batch__id=1)
        raw_path = file.get_raw_path()
        correct_path = "/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072.CEL"

        self.assertEqual(raw_path, correct_path)

    def test_temp_dir(self):
        file = File.objects.get(batch__id=1)
        temp_dir = file.get_temp_dir()

        self.assertEqual(temp_dir, "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1")

    def test_temp_pre_path(self):
        file = File.objects.get(batch__id=1)
        temp_pre_path = file.get_temp_pre_path()
        correct_path = "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1/GSM1426072.CEL"

        self.assertEqual(temp_pre_path, correct_path)

    def test_temp_post_path(self):
        file = File.objects.get(batch__id=1)
        temp_post_path = file.get_temp_post_path()
        correct_path = "/home/user/data_store/temp/A-AFFY-141/AFFY_TO_PCL/batch_1/GSM1426072.PCL"

        self.assertEqual(temp_post_path, correct_path)

    def test_processed_dir_s3(self):
        batches.USE_S3 = True
        file = File.objects.get(batch__id=1)
        processed_dir = file.get_processed_dir()

        self.assertEqual(processed_dir, "processed/A-AFFY-141/AFFY_TO_PCL")

    def test_processed_dir_no_s3(self):
        batches.USE_S3 = False
        file = File.objects.get(batch__id=1)
        processed_dir = file.get_processed_dir()

        self.assertEqual(processed_dir, "/home/user/data_store/processed/A-AFFY-141/AFFY_TO_PCL")

    def test_processed_path_s3(self):
        batches.USE_S3 = True
        file = File.objects.get(batch__id=1)
        processed_path = file.get_processed_path()
        correct_path = "processed/A-AFFY-141/AFFY_TO_PCL/GSM1426072.PCL"

        self.assertEqual(processed_path, correct_path)

    def test_processed_path_no_s3(self):
        batches.USE_S3 = False
        file = File.objects.get(batch__id=1)
        processed_path = file.get_processed_path()
        correct_path = "/home/user/data_store/processed/A-AFFY-141/AFFY_TO_PCL/GSM1426072.PCL"

        self.assertEqual(processed_path, correct_path)
