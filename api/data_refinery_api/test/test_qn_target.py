from django.urls import reverse
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_common.models import ComputationalResult, ComputationalResultAnnotation, Organism


class APITestCases(APITestCase):
    def setUp(self):
        self.homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        self.homo_sapiens.save()
        self.danio_rerio = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
        self.danio_rerio.save()

    def tearDown(self):
        Organism.objects.all().delete()

    def test_qn_endpoints(self):
        # create two qn endpoints

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

        self.assertEqual(len(response.json()), 2)
