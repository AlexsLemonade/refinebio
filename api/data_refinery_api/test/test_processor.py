from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.tests import API_VERSION
from data_refinery_common.models import (
    ComputationalResult,
    Organism,
    OrganismIndex,
    Processor,
    Sample,
    SampleResultAssociation,
)


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

    def tearDown(self):
        ComputationalResult.objects.all().delete()
        Organism.objects.all().delete()
        OrganismIndex.objects.all().delete()
        Processor.objects.all().delete()
        Sample.objects.all().delete()
        SampleResultAssociation.objects.all().delete()

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
