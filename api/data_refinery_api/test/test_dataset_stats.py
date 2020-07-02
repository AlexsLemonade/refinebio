import json

from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_common.models import (
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    Sample,
)


class APITestCases(APITestCase):
    def setUp(self):
        self.homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        self.homo_sapiens.save()

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
        sample1.is_processed = True
        sample1.organism = self.homo_sapiens
        sample1.save()

        sample2 = Sample()
        sample2.title = "2"
        sample2.accession_code = "2"
        sample2.platform_name = "ILLUMINA"
        sample2.is_processed = True
        sample2.organism = gallus_gallus
        sample2.save()

        sample3 = Sample()
        sample3.title = "3"
        sample3.accession_code = "3"
        sample3.platform_name = "ILLUMINA"
        sample3.is_processed = True
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

