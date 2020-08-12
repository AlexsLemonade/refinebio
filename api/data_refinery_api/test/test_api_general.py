import json
from unittest.mock import Mock, patch

from django.core.cache import cache
from django.core.exceptions import TooManyFieldsSent
from django.core.management import call_command
from django.http import HttpResponseForbidden, HttpResponseServerError
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.views import DatasetView, ExperimentListView
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
    OntologyTerm,
    Organism,
    OrganismIndex,
    OriginalFile,
    OriginalFileSampleAssociation,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleAttribute,
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

        # Create 26 test organisms numbered 0-25 for pagination test, so there should be 29 organisms total (with the 3 others below)
        for i in range(26):
            Organism(name=("TEST_ORGANISM_{}".format(i)), taxonomy_id=(1234 + i)).save()

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

        length = OntologyTerm()
        length.ontology_term = "PATO:0000122"
        length.human_readable_name = "length"
        length.save()

        sa = SampleAttribute()
        sa.name = length
        sa.submitter = "Refinebio Tests"
        sa.set_value(5)
        sa.sample = sample
        sa.save()

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
        experiment.update_sample_metadata_fields()
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

    def test_all_endpoints(self):
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response["X-Source-Revision"], get_env_variable("SYSTEM_VERSION"))

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

        response = self.client.get(reverse("organisms", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("organisms", kwargs={"version": API_VERSION}) + "HOMO_SAPIENS/"
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("platforms", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("institutions", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("survey_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("downloader_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        # Don't know the best way to deal with this, but since the other tests in different files
        # create objects which are then deleted, the new objects from these tests will have different
        # IDs. In this case, since this file is ran first, the IDs are 1, but this may be a problem
        # in the future.
        response = self.client.get(
            reverse("downloader_jobs", kwargs={"version": API_VERSION}) + "1/"  # change back
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("processor_jobs", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("processor_jobs", kwargs={"version": API_VERSION}) + "1/"
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("results", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("results", kwargs={"version": API_VERSION}) + "1/")
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("schema_redoc", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("transcriptome_indices", kwargs={"version": API_VERSION})
            + "?organism__name=DANIO_RERIO"
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("transcriptome_indices", kwargs={"version": API_VERSION}) + "?result_id=1"
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("search", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(
            reverse("transcriptome_indices", kwargs={"version": API_VERSION})
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse("create_dataset", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)

    def test_experiment_multiple_accessions(self):
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION})
            + "?accession_code=GSE000&accession_code=GSE123",
            follow=True,
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 2)

    # Test the query the front-end uses to find the experiment with a given
    # accession or alternate accession
    def test_experiment_alternate_accession(self):
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION})
            + "?search=alternate_accession_code:E-GEOD-000"
            + "?search=accession_code:E-GEOD-000",
            follow=True,
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 1)
        self.assertEqual(response.json()["results"][0]["alternate_accession_code"], "E-GEOD-000")

    def test_experiment_external_metadata(self):
        """Test if we can filter based on metadata supplied by an external contributor"""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}) + "?sample_metadata_fields=length",
            follow=True,
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.json()["results"]), 1)
        self.assertEqual(response.json()["results"][0]["accession_code"], "GSE123")

    def test_sample_multiple_accessions(self):
        response = self.client.get(
            reverse("samples", kwargs={"version": API_VERSION}) + "?accession_codes=123,789",
            follow=True,
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 2)

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

    def test_organism_pagination(self):
        response = self.client.get(reverse("organisms", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 25)

        # First organism on second page should be TEST_ORGANISM_25, and since 29 organisms have been created, there should be 4 on the 2nd page
        response = self.client.get(
            reverse("organisms", kwargs={"version": API_VERSION}), {"offset": 25}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()["results"]), 4)
        self.assertEqual(response.json()["results"][0]["name"], "TEST_ORGANISM_25")

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

    def test_sample_detail_experiment_accessions(self):
        response = self.client.get(
            reverse("samples_detail", kwargs={"version": API_VERSION, "accession_code": "789"})
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()["experiment_accession_codes"], ["GSE123"])

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

    def test_create_token(self):
        # First, try activating right away
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}),
            json.dumps({"is_activated": True}),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 201)
        token_id = response.json()["id"]

        # Now activate using a second request
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}), content_type="application/json"
        )
        self.assertEqual(response.status_code, 201)
        token = response.json()
        token["is_activated"] = True
        token_id = token["id"]
        response = self.client.put(
            reverse("token_id", kwargs={"id": token_id, "version": API_VERSION}),
            json.dumps(token),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)

        activated_token = response.json()
        self.assertEqual(activated_token["id"], token_id)
        self.assertEqual(activated_token["is_activated"], True)
