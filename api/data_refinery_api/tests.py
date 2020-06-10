import json
from unittest.mock import Mock, patch

from django.core.cache import cache
from django.core.management import call_command
from django.http import HttpResponseForbidden, HttpResponseServerError
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.views import DatasetView, ExperimentList
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentOrganismAssociation,
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
from data_refinery_common.models.documents import ExperimentDocument
from data_refinery_common.utils import get_env_variable

API_VERSION = "v1"


class APITestCases(APITestCase):
    def setUp(self):
        # Saving this for if we have protected endpoints
        # self.superuser = User.objects.create_superuser('john', 'john@snow.com', 'johnpassword')
        # self.client.login(username='john', password='johnpassword')
        # self.user = User.objects.create(username="mike")

        experiment = Experiment()
        experiment.accession_code = "GSE000"
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
        """ Good bye """
        Experiment.objects.all().delete()
        ExperimentAnnotation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()

    def test_all_endpoints(self):
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response["X-Source-Revision"], get_env_variable("SYSTEM_VERSION"))

        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("samples", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}),
            {"ids": str(self.sample.id) + ",1000"},
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}),
            {"accession_codes": str(self.sample.accession_code) + ",1000"},
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("samples", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("organisms", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("platforms", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("institutions", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("survey_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("downloader_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("processor_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("results", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("results", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("schema_redoc", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("search", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("transcriptome_indices", kwargs={"version": API_VERSION})
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("create_dataset", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)

    def test_sample_pagination(self):
        response = self.client.get(reverse("samples", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 2)

        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}), {"limit": 1}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 1)

        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}), {"limit": 1, "ordering": "-title"}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()["results"][0]["title"], "789")

        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}), {"limit": 1, "ordering": "title"}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()["results"][0]["title"], "123")

    def test_fetching_experiment_samples(self):
        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}),
            {"experiment_accession_code": self.experiment.accession_code},
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 1)
        self.assertEqual(response.json()["results"][0]["accession_code"], "789")

        # Expect 404 if the experiment accession code isn't valid
        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}),
            {"experiment_accession_code": "wrong-accession-code"},
        )
        self.assertEqual(response.status_code, 404)

    def test_fetching_organism_index(self):
        organism_index_id = OrganismIndex.objects.all().first().id
        response = self.client.get(
            reverse(
                "transcriptome_indices_read",
                kwargs={"id": organism_index_id, "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()["index_type"], "TRANSCRIPTOME_SHORT")

        # Expect 404 if the transcriptome index id is not valid
        response = self.client.get(
            reverse("transcriptome_indices_read", kwargs={"id": 0, "version": API_VERSION})
        )
        self.assertEqual(response.status_code, 404)

    def test_compendia(self):

        result = ComputationalResult()
        result.save()

        hsc1 = ComputedFile()
        hsc1.absolute_file_path = "/null/1.tsv"
        hsc1.filename = "1.tsv"
        hsc1.sha1 = "abc"
        hsc1.size_in_bytes = 1
        hsc1.is_smashable = False
        hsc1.is_qn_target = False
        hsc1.result = result
        hsc1.is_compendia = True
        hsc1.compendia_organism = self.homo_sapiens
        hsc1.compendia_version = 1
        hsc1.s3_bucket = "dr-compendia"
        hsc1.s3_key = "hsc1.tsv"
        hsc1.save()

        hsc2 = ComputedFile()
        hsc2.absolute_file_path = "/null/2.tsv"
        hsc2.filename = "2.tsv"
        hsc2.sha1 = "abc"
        hsc2.size_in_bytes = 1
        hsc2.is_smashable = False
        hsc2.is_qn_target = False
        hsc2.result = result
        hsc2.is_compendia = True
        hsc2.compendia_organism = self.homo_sapiens
        hsc2.compendia_version = 2
        hsc2.s3_bucket = "dr-compendia"
        hsc2.s3_key = "hsc2.tsv"
        hsc2.save()

        drc1 = ComputedFile()
        drc1.absolute_file_path = "/null/1.tsv"
        drc1.filename = "1.tsv"
        drc1.sha1 = "abc"
        drc1.size_in_bytes = 1
        drc1.is_smashable = False
        drc1.is_qn_target = False
        drc1.result = result
        drc1.is_compendia = True
        drc1.compendia_organism = self.danio_rerio
        drc1.compendia_version = 1
        drc1.s3_bucket = "dr-compendia"
        drc1.s3_key = "drc2.tsv"
        drc1.save()

        response = self.client.get(
            reverse("computed_files", kwargs={"version": API_VERSION}), {"is_compendia": True}
        )
        response_json = response.json()["results"]
        self.assertEqual(3, len(response_json))
        # Prove that the download_url field is missing and not None.
        self.assertEqual("NotPresent", response_json[0].get("download_url", "NotPresent"))

        # We don't actually want AWS to generate a temporary URL for
        # us, and it won't unless we're running in the cloud, but if
        # we provide an API Token and use the WithUrl serializer then
        # it will set the download_url field to None rather than
        # generate one.

        # Create a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}), content_type="application/json"
        )
        token = response.json()
        token["is_activated"] = True
        token_id = token["id"]
        response = self.client.put(
            reverse("token_id", kwargs={"id": token_id, "version": API_VERSION}),
            json.dumps(token),
            content_type="application/json",
        )

        response = self.client.get(
            reverse("computed_files", kwargs={"version": API_VERSION}),
            {"is_compendia": True},
            HTTP_API_KEY=token_id,
        )
        response_json = response.json()["results"]
        self.assertEqual(3, len(response_json))
        self.assertIsNone(response_json[0]["download_url"])

    @patch("data_refinery_common.message_queue.send_job")
    def test_create_dataset_fails_without_email(self, mock_send_job):

        # Get a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}), content_type="application/json"
        )
        token = response.json()
        token["is_activated"] = True
        token_id = token["id"]
        response = self.client.put(
            reverse("token_id", kwargs={"id": token_id, "version": API_VERSION}),
            json.dumps(token),
            content_type="application/json",
        )

        activated_token = response.json()
        self.assertEqual(activated_token["id"], token_id)
        self.assertEqual(activated_token["is_activated"], True)

        # Good, except for missing email.
        jdata = json.dumps({"data": {"A": ["B"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
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

    @patch("data_refinery_common.message_queue.send_job")
    def test_create_update_dataset(self, mock_send_job):

        # Get a token first
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}), content_type="application/json"
        )
        token = response.json()
        token["is_activated"] = True
        token_id = token["id"]
        response = self.client.put(
            reverse("token_id", kwargs={"id": token_id, "version": API_VERSION}),
            json.dumps(token),
            content_type="application/json",
        )

        activated_token = response.json()
        self.assertEqual(activated_token["id"], token_id)
        self.assertEqual(activated_token["is_activated"], True)
        # Good
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"A": ["B"]}})
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
        self.assertEqual(response.json()["data"]["A"], ["B"])

        # Bad (Duplicates)
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"A": ["B", "B", "B"]}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Bad (Empty Experiment)
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"A": []}})
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)

        # Bad (Empty Experiment)
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {}})
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
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"A": ["C"]}})
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["id"], good_id)
        self.assertEqual(response.json()["data"], json.loads(jdata)["data"])
        self.assertEqual(response.json()["data"]["A"], ["C"])

        # Can't update if started
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = True
        dataset.save()
        jdata = json.dumps({"email_address": "baz@gmail.com", "data": {"A": ["D"]}})
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertNotEqual(response.json()["data"]["A"], ["D"])

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
                "data": {"A": ["D"]},
                "start": True,
                "no_send_job": True,
                "token_id": "HEYO",
            }
        )
        response = self.client.put(
            reverse("dataset", kwargs={"id": good_id, "version": API_VERSION}),
            jdata,
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 500)

        jdata = json.dumps(
            {
                "data": {"A": ["D"]},
                "start": True,
                "no_send_job": True,
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
                "data": {"A": ["D"]},
                "start": True,
                "no_send_job": True,
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

    def test_processed_samples_only(self):
        """ Don't return unprocessed samples """
        experiment = Experiment()
        experiment.accession_code = "GSX12345"
        experiment.is_public = True
        experiment.save()

        sample = Sample()
        sample.title = "I am unprocessed"
        sample.accession_code = "GSXUnprocessed"
        sample.is_processed = False
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        # we return all experiments
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}), {"search": "GSX12345"}
        )
        self.assertEqual(response.json()["count"], 1)

        # check requesting only experiments with processed samples
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "GSX12345", "num_processed_samples__gt": 0},
        )
        self.assertEqual(response.json()["count"], 0)

        sample2 = Sample()
        sample2.title = "I am processed"
        sample2.accession_code = "GSXProcessed"
        sample2.is_processed = True
        sample2.save()

        experiment_sample2_association = ExperimentSampleAssociation()
        experiment_sample2_association.sample = sample2
        experiment_sample2_association.experiment = experiment
        experiment_sample2_association.save()

        # update cached values
        experiment.num_total_samples = 2
        experiment.num_processed_samples = 1
        experiment.save()

        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}), {"search": "GSX12345"}
        )
        self.assertEqual(response.json()["count"], 1)

        self.assertEqual(len(experiment.processed_samples), 1)

        experiment.delete()
        sample.delete()
        sample2.delete()

    def test_dataset_stats(self):
        """ Test the dataset stats endpoint """

        gallus_gallus = Organism(name="GALLUS_GALLUS", taxonomy_id=9031, is_scientific_name=True)
        gallus_gallus.save()
        equus_ferus = Organism(name="EQUUS_FERUS", taxonomy_id=1114792, is_scientific_name=True)
        equus_ferus.save()

        ex = Experiment()
        ex.accession_code = "XYZ123"
        ex.title = "XYZ123"
        ex.description = "XYZ123"
        ex.technology = "MICROARRAY"
        ex.submitter_institution = "XYZ123"
        ex.save()

        ex2 = Experiment()
        ex2.accession_code = "ABC789"
        ex2.title = "ABC789"
        ex2.description = "ABC789"
        ex2.technology = "RNA-SEQ"
        ex2.submitter_institution = "Funkytown"
        ex2.save()

        sample1 = Sample()
        sample1.title = "1"
        sample1.accession_code = "1"
        sample1.platform_name = "AFFY"
        sample1.organism = self.homo_sapiens
        sample1.save()

        sample2 = Sample()
        sample2.title = "2"
        sample2.accession_code = "2"
        sample2.platform_name = "ILLUMINA"
        sample2.organism = gallus_gallus
        sample2.save()

        sample3 = Sample()
        sample3.title = "3"
        sample3.accession_code = "3"
        sample3.platform_name = "ILLUMINA"
        sample3.organism = gallus_gallus
        sample3.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment = ex
        xoa.organism = self.homo_sapiens
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment = ex2
        xoa.organism = gallus_gallus
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment = ex2
        xoa.organism = equus_ferus
        xoa.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample1
        experiment_sample_association.experiment = ex
        experiment_sample_association.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample2
        experiment_sample_association.experiment = ex2
        experiment_sample_association.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample3
        experiment_sample_association.experiment = ex2
        experiment_sample_association.save()

        jdata = json.dumps(
            {"email_address": "baz@gmail.com", "data": {"XYZ123": ["1"], "ABC789": ["2"]}}
        )
        response = self.client.post(
            reverse("create_dataset", kwargs={"version": API_VERSION}),
            jdata,
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        self.assertEqual(response.json()["data"], json.loads(jdata)["data"])
        good_id = response.json()["id"]

        # Check that we can fetch these sample details via samples API
        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}), {"dataset_id": good_id}
        )
        self.assertEqual(response.json()["count"], 2)

    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_ok(self, mock_client):
        # We don't even import raven if it's a good response.
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        mock_client.is_enabled.assert_not_called()

    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_404(self, mock_client):
        # We don't send anything to raven if it's not enabled
        mock_client.is_enabled.side_effect = lambda: False
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )[:-1]
            + "aasdas/"
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

    @patch.object(ExperimentList, "get")
    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_403(self, mock_client, mock_get_method):
        mock_get_method.side_effect = Mock(return_value=HttpResponseForbidden())
        # A 403 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = Mock(return_value=True)
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 403)
        mock_client.captureMessage.assert_called()

    @patch.object(ExperimentList, "get")
    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_500(self, mock_client, mock_get_method):
        def raise_error(_):
            raise KeyError()

        mock_get_method.side_effect = Mock(return_value=HttpResponseServerError())
        # A 500 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = Mock(return_value=True)
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 500)
        mock_client.captureMessage.assert_called()

    def test_qn_endpoints(self):

        # create two additional qn endpoints

        result = ComputationalResult()
        result.commands.append("create_qn_target.py")
        result.is_ccdl = True
        result.is_public = True
        result.processor = None
        result.save()

        cra = ComputationalResultAnnotation()
        cra.result = result
        cra.data = {
            "organism_id": self.danio_rerio.id,  # Danio
            "is_qn": True,
            "platform_accession_code": "zebrafish",
            "samples": [],
            "geneset": str(["RWWJ000001", "RWWJ000002"]),
        }
        cra.save()
        cra = ComputationalResultAnnotation()
        cra.result = result
        cra.data = {
            "organism_id": self.homo_sapiens.id,  # IDK
            "is_qn": True,
            "platform_accession_code": "zebrafishplusone",
            "samples": [],
            "geneset": str(["RWWJ000003", "RWWJ000004"]),
        }
        cra.save()

        self.homo_sapiens.qn_target = result
        self.homo_sapiens.save()
        self.danio_rerio.qn_target = result
        self.danio_rerio.save()

        response = self.client.get(reverse("qn_targets_available", kwargs={"version": API_VERSION}))
        # there's another qn endpoint that is created in the setup method of this test case
        self.assertEqual(len(response.json()), 3)


MOCK_NOMAD_RESPONSE = [
    {
        "CreateIndex": 5145,
        "ID": "TXIMPORT",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 0},
            "CreateIndex": 5145,
            "JobID": "TXIMPORT",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "TXIMPORT",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_1_2323",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 1},
            "CreateIndex": 5145,
            "JobID": "SALMON_1_2323",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_1_2323",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_3_2534/dispatch-213123-123123-123123",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 0},
            "CreateIndex": 5145,
            "JobID": "SALMON_3_2534/dispatch-213123-123123-123123",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_3_2534/dispatch-213123-123123-123123",
        "ParameterizedJob": True,
        "ParentID": "SALMON_1_2323",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_2_2323",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 1},
            "CreateIndex": 5145,
            "JobID": "SALMON_1_2323",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_1_2323",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
]


class StatsTestCases(APITestCase):
    @patch("data_refinery_common.utils.get_nomad_jobs", return_value=[])
    def test_nomad_stats_empty(self, mock_get_nomad_jobs):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["nomad_running_jobs"], 0)
        self.assertEqual(response.json()["nomad_pending_jobs"], 0)

    @patch("data_refinery_common.utils.get_nomad_jobs", return_value=MOCK_NOMAD_RESPONSE)
    def test_nomad_stats(self, mock_get_nomad_jobs):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["nomad_running_jobs"], 2)
        self.assertEqual(response.json()["nomad_pending_jobs"], 0)
        self.assertEqual(response.json()["nomad_running_jobs_by_type"]["SALMON"], 2)
        self.assertEqual(response.json()["nomad_running_jobs_by_volume"]["1"], 1)
        self.assertEqual(response.json()["nomad_running_jobs_by_volume"]["2"], 1)


class ESTestCases(APITestCase):
    @classmethod
    def setUpClass(cls):
        super(ESTestCases, cls).setUpClass()  # ref https://stackoverflow.com/a/29655301/763705

        """Set up class."""
        experiment = Experiment()
        experiment.accession_code = "GSE000-X"
        experiment.title = "NONONONO"
        experiment.description = "Boooooourns. Wasabi."
        experiment.technology = "RNA-SEQ"
        experiment.save()

        experiment = Experiment()
        experiment.accession_code = "GSE123-X"
        experiment.title = "Hey Ho Let's Go"
        experiment.description = (
            "This is a very exciting test experiment. Faygo soda. Blah blah blah."
        )
        experiment.technology = "MICROARRAY"
        experiment.num_processed_samples = 1  # added below
        experiment.num_total_samples = 1
        experiment.num_downloadable_samples = 1
        experiment.save()

        experiment_annotation = ExperimentAnnotation()
        experiment_annotation.data = {"hello": "world", "123": 456}
        experiment_annotation.experiment = experiment
        experiment_annotation.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.save()

        organism = Organism(
            name="AILUROPODA_MELANOLEUCA", taxonomy_id=9646, is_scientific_name=True
        )
        organism.save()

        sample = Sample()
        sample.title = "789"
        sample.accession_code = "789"
        sample.is_processed = True
        sample.organism = organism
        sample.save()

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

        # associate the experiment with the sample
        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        result = ComputationalResult()
        result.save()

        # and create a qn tarjet for the sample
        computational_result = ComputationalResultAnnotation()
        computational_result.result = result
        computational_result.data = {"is_qn": True, "organism_id": sample.organism.id}
        computational_result.save()

        # and associate it with the sample organism
        sample.organism.qn_target = result
        sample.organism.save()

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

        # clear default cache and reindex
        # otherwise the organisms with qn_targes will be cached.
        cache.clear()
        call_command("search_index", "--rebuild", "-f")

    def test_es_endpoint(self):
        """ Test basic ES functionality  """
        response = self.client.get(reverse("search", kwargs={"version": API_VERSION}))

        """ Test basic ES functionality """
        es_search_result = ExperimentDocument.search().filter("term", description="soda")
        es_search_result_qs = es_search_result.to_queryset()
        self.assertEqual(len(es_search_result_qs), 1)

        # Sanity
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 2)

        # test sample counts in filters
        self.assertEqual(response.json()["facets"]["has_publication"]["false"], 1)
        self.assertEqual(response.json()["facets"]["technology"]["microarray"], 1)
        self.assertEqual(response.json()["facets"]["technology"]["rna-seq"], 0)

        # Basic Search
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}), {"search": "soda"}
        )
        self.assertEqual(response.json()["count"], 1)

        # Positive filter result
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "soda", "technology": "microarray"},
        )
        self.assertEqual(response.json()["count"], 1)

        # Negative filter result
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "soda", "technology": "rna"},
        )
        self.assertEqual(response.json()["count"], 0)


class ProcessorTestCases(APITestCase):
    def setUp(self):
        salmon_quant_env = {
            "os_distribution": "Ubuntu 16.04.4 LTS",
            "os_pkg": {"python3": "3.5.1-3", "python3-pip": "8.1.1-2ubuntu0.4"},
            "cmd_line": {"salmon --version": "salmon 0.9.1"},
            "python": {"Django": "2.0.6", "data-refinery-common": "0.5.0"},
        }
        self.salmon_quant_proc = Processor.objects.create(
            name="Salmon Quant",
            version="0.45",
            docker_image="ccdl/salmon_img:v1.23",
            environment=salmon_quant_env,
        )

        salmontools_env = {
            "os_distribution": "Ubuntu 16.04.4 LTS",
            "os_pkg": {
                "python3": "3.5.1-3",
                "python3-pip": "8.1.1-2ubuntu0.4",
                "g++": "4:5.3.1-1ubuntu1",
                "cmake": "3.5.1-1ubuntu3",
            },
            "cmd_line": {"salmontools --version": "Salmon Tools 0.1.0"},
            "python": {"Django": "2.0.6", "data-refinery-common": "0.5.0"},
        }
        Processor.objects.create(
            name="Salmontools",
            version="1.83",
            docker_image="ccdl/salmontools_img:v0.45",
            environment=salmontools_env,
        )

    def test_endpoint(self):
        response = self.client.get(reverse("processors", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        processors = response.json()["results"]

        self.assertEqual(processors[0]["name"], "Salmon Quant")
        self.assertEqual(processors[0]["environment"]["os_pkg"]["python3"], "3.5.1-3")

        self.assertEqual(processors[1]["name"], "Salmontools")
        self.assertEqual(
            processors[1]["environment"]["cmd_line"]["salmontools --version"], "Salmon Tools 0.1.0"
        )

    def test_processor_and_organism_in_sample(self):
        sample = Sample.objects.create(accession_code="ACCESSION", title="fake sample")
        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        homo_sapiens.save()
        transcriptome_result = ComputationalResult.objects.create()
        organism_index = OrganismIndex.objects.create(
            organism=homo_sapiens, result=transcriptome_result, index_type="TRANSCRIPTOME_LONG"
        )
        result = ComputationalResult.objects.create(
            processor=self.salmon_quant_proc, organism_index=organism_index
        )
        SampleResultAssociation.objects.create(sample=sample, result=result)

        response = self.client.get(
            reverse(
                "samples_detail",
                kwargs={"accession_code": sample.accession_code, "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        processor = response.json()["results"][0]["processor"]
        self.assertEqual(processor["name"], self.salmon_quant_proc.name)
        self.assertEqual(
            processor["environment"]["os_pkg"]["python3"],
            self.salmon_quant_proc.environment["os_pkg"]["python3"],
        )

        organism_index = response.json()["results"][0]["organism_index"]
        self.assertEqual(organism_index["result_id"], transcriptome_result.id)
        self.assertEqual(organism_index["index_type"], "TRANSCRIPTOME_LONG")
