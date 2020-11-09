import json
from unittest.mock import Mock, patch

from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    OriginalFile,
    OriginalFileSampleAssociation,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleResultAssociation,
)


class DatasetTestCase(APITestCase):
    def setUp(self):
        # Saving this for if we have protected endpoints
        # self.superuser = User.objects.create_superuser('john', 'john@snow.com', 'johnpassword')
        # self.client.login(username='john', password='johnpassword')
        # self.user = User.objects.create(username="mike")

        experiment = Experiment()
        experiment.accession_code = "GSE000"
        experiment.alternate_accession_code = "E-GEOD-000"
        experiment.title = "NONONONO"
        experiment.description = "Boooooourns. Wasabi."
        experiment.technology = "RNA-SEQ"
        experiment.save()

        experiment = Experiment()
        experiment.accession_code = "GSE123"
        experiment.title = "Hey Ho Let's Go"
        experiment.description = (
            "This is a very exciting test experiment. Faygo soda. Blah blah blah."
        )
        experiment.technology = "MICROARRAY"
        experiment.save()
        self.experiment = experiment

        experiment_annotation = ExperimentAnnotation()
        experiment_annotation.data = {"hello": "world", "123": 456}
        experiment_annotation.experiment = experiment
        experiment_annotation.save()

        ailuropoda = Organism(
            name="AILUROPODA_MELANOLEUCA", taxonomy_id=9646, is_scientific_name=True
        )
        ailuropoda.save()
        self.homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        self.homo_sapiens.save()
        self.danio_rerio = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
        self.danio_rerio.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.is_processed = True
        sample.organism = ailuropoda
        sample.save()

        sample = Sample()
        sample.title = "789"
        sample.accession_code = "789"
        sample.is_processed = True
        sample.organism = ailuropoda
        sample.save()
        self.sample = sample

        # add qn target for sample organism
        result = ComputationalResult()
        result.commands.append("create_qn_target.py")
        result.is_ccdl = True
        result.is_public = True
        result.processor = None
        result.save()

        cra = ComputationalResultAnnotation()
        cra.result = result
        cra.data = {"organism_id": ailuropoda.id, "is_qn": True}
        cra.save()

        ailuropoda.qn_target = result
        ailuropoda.save()

        sample_annotation = SampleAnnotation()
        sample_annotation.data = {"goodbye": "world", "789": 123}
        sample_annotation.sample = sample
        sample_annotation.save()

        original_file = OriginalFile()
        original_file.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.sample = sample
        original_file_sample_association.original_file = original_file
        original_file_sample_association.save()

        downloader_job = DownloaderJob()
        downloader_job.save()

        download_assoc = DownloaderJobOriginalFileAssociation()
        download_assoc.original_file = original_file
        download_assoc.downloader_job = downloader_job
        download_assoc.save()

        processor_job = ProcessorJob()
        processor_job.save()

        processor_assoc = ProcessorJobOriginalFileAssociation()
        processor_assoc.original_file = original_file
        processor_assoc.processor_job = processor_job
        processor_assoc.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()
        experiment.num_total_samples = 1
        experiment.num_processed_samples = 1
        experiment.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        processor = Processor()
        processor.name = "Salmon Quant"
        processor.version = "v9.9.9"
        processor.docker_image = "dr_salmon"
        processor.environment = '{"some": "environment"}'
        processor.save()

        computational_result_short = ComputationalResult(processor=processor)
        computational_result_short.save()

        organism_index = OrganismIndex()
        organism_index.index_type = "TRANSCRIPTOME_SHORT"
        organism_index.organism = self.danio_rerio
        organism_index.result = computational_result_short
        organism_index.absolute_directory_path = (
            "/home/user/data_store/salmon_tests/TRANSCRIPTOME_INDEX/SHORT"
        )
        organism_index.is_public = True
        organism_index.s3_url = "not_blank"
        organism_index.save()

        return

    def tearDown(self):
        Experiment.objects.all().delete()
        ExperimentSampleAssociation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()

    @patch("data_refinery_api.views.dataset.send_job", lambda *args: True)
    def test_create_dataset_fails_without_email(self):
        # Create a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}),
            json.dumps({"is_activated": True}),
            content_type="application/json",
        )
        token_id = response.json()["id"]

        # Good, except for missing email.
        jdata = json.dumps({"start": True, "data": {"GSE123": ["789"]}, "token_id": token_id})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Good, except for empty email.
        jdata = json.dumps(
            {"start": True, "data": {"GSE123": ["789"]}, "token_id": token_id, "email_address": "",}
        )
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # You should not have to provide an email until you set start=True
        jdata = json.dumps({"data": {"GSE123": ["789"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)

        # Good, except for missing email.
        jdata = json.dumps({"start": True, "data": {"GSE123": ["789"]}, "token_id": token_id})
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Good, except for invalid email.
        jdata = json.dumps({"email_address": "bad format!", "data": {"A": ["B"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

    @patch("data_refinery_api.views.dataset.send_job", lambda *args: True)
    def test_starting_dataset_fails_without_experiments(self):

        # Create a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}),
            json.dumps({"is_activated": True}),
            content_type="application/json",
        )

        token_id = response.json()["id"]

        # Good, except for having zero experiments
        data = {"email_address": "baz@gmail.com", "data": {}}
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            json.dumps(data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)

        data["start"] = True
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps({"start": True}),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)

        # On the other hand, if we add experiments before we start processing, this should be okay
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)

        data = {
            **data,
            "start": True,
            "data": {"GSE123": ["789"]},
            "token_id": token_id,
            "no_send_job": True,
        }
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps(data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)

    def test_dataset_adding_non_downloadable_samples_fails(self):
        # Make a sample that is not downloadable
        sample1 = Sample()
        sample1.title = "456"
        sample1.accession_code = "456"
        sample1.platform_name = "AFFY"
        sample1.is_processed = False
        sample1.organism = self.homo_sapiens
        sample1.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample1
        experiment_sample_association.experiment = self.experiment
        experiment_sample_association.save()

        # Bad, 456 is not processed
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["456"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(
            "Non-downloadable sample(s) in dataset", response.json()["message"],
        )
        self.assertEqual(response.json()["details"], ["456"])

        # Bad, 567 does not exist
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["567"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertIn(
            "Sample(s) in dataset do not exist on refine", response.json()["message"],
        )
        self.assertEqual(response.status_code, 400)

        # Good, 789 is processed
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["789"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)

        # Bad, 456 does not have a quant.sf file
        post_data = {"email_address": "baz@gmail.com", "data": {}}
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            json.dumps(post_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)

        put_data = {**post_data, "data": {"GSE123": ["456"]}, "quant_sf_only": True}
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps(put_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(
            "Sample(s) in dataset are missing quant.sf files", response.json()["message"],
        )
        self.assertEqual(response.json()["details"], ["456"])

        # Bad, none of the samples in GSE123 have a quant.sf file
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            json.dumps(post_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps({**put_data, "data": {"GSE123": ["ALL"]}}),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn(
            "Experiment(s) in dataset have zero downloadable samples", response.json()["message"],
        )
        self.assertEqual(response.json()["details"], ["GSE123"])

        # Make 456 have a quant.sf file
        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample1
        sra.result = result
        sra.save()

        computed_file = ComputedFile()
        computed_file.s3_key = "smasher-test-quant.sf"
        computed_file.s3_bucket = "data-refinery-test-assets"
        computed_file.filename = "quant.sf"
        computed_file.result = result
        computed_file.size_in_bytes = 42
        computed_file.save()

        # Good, 456 does have a quant.sf file
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            json.dumps(post_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)

        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps(put_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)

        # Good, a sample in GSE123 has a quant.sf file
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            json.dumps(post_data),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)
        response = self.client.put(
            reverse("dataset", kwargs={"id": response.json()["id"], "version": API_VERSION}),
            json.dumps({**put_data, "data": {"GSE123": ["ALL"]}}),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)

    @patch("data_refinery_api.views.dataset.send_job", lambda *args: True)
    def test_create_update_dataset(self):

        # Create a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}),
            json.dumps({"is_activated": True}),
            content_type="application/json",
        )
        token_id = response.json()["id"]

        # Good
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["789"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        self.assertEqual(response.json()["data"], json.loads(jdata)["data"])
        good_id = response.json()["id"]

        response = self.client.get(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION})
        )
        self.assertEqual(response.json()["id"], good_id)
        self.assertEqual(response.json()["data"], json.loads(jdata)["data"])
        self.assertEqual(response.json()["data"]["GSE123"], ["789"])

        # Bad (Duplicates)
        jdata = json.dumps(
            {"email_address": "baz@gmail.com", "data": {"GSE123": ["789", "789", "789"]}}
        )
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Bad (Empty Experiment)
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": []}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Update, just an experiment accession
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["ALL"]}})
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["id"], good_id)
        # We are not entirely RESTful here, that's okay.
        self.assertNotEqual(response.json()["data"], json.loads(jdata)["data"])
        self.assertEqual(response.json()["data"]["GSE123"], ["789"])

        # Update
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["789"]}})
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["id"], good_id)
        self.assertEqual(response.json()["data"], json.loads(jdata)["data"])
        self.assertEqual(response.json()["data"]["GSE123"], ["789"])

        # Can't update if started
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = True
        dataset.save()
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"GSE123": ["123"]}})
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        # The request should succeed, but the dataset should not be modified
        self.assertEqual(response.status_code, 200)
        self.assertNotEqual(response.json()["data"]["GSE123"], ["123"])

        # Bad
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": 123})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)

        # This will actually kick off a job if we don't patch send_job or supply no_send_job
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = False
        dataset.save()

        # With bad token first
        jdata = json.dumps(
            {
                "email_address": "baz@gmail.com",
                "data": {"GSE123": ["789"]},
                "start": True,
                "token_id": "HEYO",
            }
        )
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 400)

        jdata = json.dumps(
            {
                "data": {"GSE123": ["789"]},
                "start": True,
                "token_id": token_id,
                "email_address": "trust@verify.com",
                "email_ccdl_ok": True,
            }
        )
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.json()["is_processing"], True)

        ds = Dataset.objects.get(id=response.json()["id"])
        self.assertEqual(ds.email_address, "trust@verify.com")
        self.assertTrue(ds.email_ccdl_ok)

        # Reset the dataset so we can test providing an API token via
        # HTTP Header.
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = False
        dataset.email_address = "baz@gmail.com"
        dataset.email_ccdl_ok = False
        dataset.save()

        jdata = json.dumps(
            {
                "data": {"GSE123": ["789"]},
                "start": True,
                "email_address": "trust@verify.com",
                "email_ccdl_ok": True,
            }
        )
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
            HTTP_API_KEY=token_id,
        )

        self.assertEqual(response.json()["is_processing"], True)

        ds = Dataset.objects.get(id=response.json()["id"])
        self.assertEqual(ds.email_address, "trust@verify.com")
        self.assertTrue(ds.email_ccdl_ok)

        # Test creating and starting a dataset in the same action
        jdata = json.dumps(
            {
                "data": {"GSE123": ["789"]},
                "start": True,
                "email_address": "trust@verify.com",
                "email_ccdl_ok": True,
            }
        )
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
            HTTP_API_KEY=token_id,
        )
        self.assertEqual(response.json()["is_processing"], True)

        ds = Dataset.objects.get(id=response.json()["id"])
        self.assertEqual(ds.email_address, "trust@verify.com")
        self.assertTrue(ds.email_ccdl_ok)
