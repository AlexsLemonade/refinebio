import json

from django.urls import reverse
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_common.models import ComputationalResult, ComputedFile, Organism


class APITestCases(APITestCase):
    def setUp(self):
        self.homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        self.homo_sapiens.save()
        self.danio_rerio = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
        self.danio_rerio.save()
        self.rattus_rattus = Organism(
            name="RATTUS_RATTUS", taxonomy_id=10117, is_scientific_name=True
        )
        self.rattus_rattus.save()

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
        hsc1.quant_sf_only = False
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
        hsc2.quant_sf_only = False
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
        drc1.quant_sf_only = True
        drc1.compendia_organism = self.danio_rerio
        drc1.compendia_version = 1
        drc1.s3_bucket = "dr-compendia"
        drc1.s3_key = "drc2.tsv"
        drc1.save()

        rrc1 = ComputedFile()
        rrc1.absolute_file_path = "/null/1.tsv"
        rrc1.filename = "1.tsv"
        rrc1.sha1 = "abc"
        rrc1.size_in_bytes = 1
        rrc1.is_smashable = False
        rrc1.is_qn_target = False
        rrc1.result = result
        rrc1.is_compendia = False
        rrc1.quant_sf_only = False
        rrc1.compendia_organism = self.rattus_rattus
        rrc1.compendia_version = 1
        rrc1.s3_bucket = "dr-compendia"
        rrc1.s3_key = "drc2.tsv"
        rrc1.save()

    def tearDown(self):
        ComputedFile.objects.all().delete()

    def test_compendia(self):
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

        # create token
        response = self.client.post(
            reverse("token", kwargs={"version": API_VERSION}),
            json.dumps({"is_activated": True}),
            content_type="application/json",
        )
        token_id = response.json()["id"]

        response = self.client.get(
            reverse("computed_files", kwargs={"version": API_VERSION}),
            {"is_compendia": True},
            HTTP_API_KEY=token_id,
        )
        response_json = response.json()["results"]
        self.assertEqual(3, len(response_json))
        self.assertIsNone(response_json[0]["download_url"])

    def test_organism_has_compendia(self):
        response = self.client.get(
            reverse("organisms", kwargs={"version": API_VERSION}), {"has_compendia": True},
        )
        response_json = response.json()["results"]
        self.assertEqual(2, len(response_json))

        response = self.client.get(
            reverse("organisms", kwargs={"version": API_VERSION}), {"has_compendia": False},
        )
        response_json = response.json()["results"]
        self.assertEqual(1, len(response_json))
