import json

from django.urls import reverse
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Organism,
)


class APITestCases(APITestCase):
    def setUp(self):
        self.homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        self.homo_sapiens.save()
        self.danio_rerio = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
        self.danio_rerio.save()
    
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
    

