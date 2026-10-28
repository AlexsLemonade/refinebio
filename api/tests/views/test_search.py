import json

from django.core.cache import cache
from django.core.exceptions import TooManyFieldsSent
from django.core.management import call_command
from django.urls import reverse
from rest_framework.test import APITestCase

from tests.views.test_api_general import API_VERSION

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleResultAssociation,
)
from data_refinery_common.models.documents import ExperimentDocument

MOUSE = "Mus musculus"

ECOLI_STRAIN_NAME = "Escherichia coli str. k-12 substr. mg1655"


class ESTestCases(APITestCase):
    # Fixture codes shared across the extended tests below.
    # EXP_A — microarray, has_publication=True,  10 downloadable, human
    # EXP_B — rna-seq,   has_publication=False,  5 downloadable, human
    # EXP_C — rna-seq,   has_publication=True,  20 downloadable, mouse
    EXP_A = "GSEPA1"
    EXP_B = "GSEPA2"
    EXP_C = "GSEPA3"

    @classmethod
    def _make_qn_organism(cls, name, taxonomy_id):
        organism = Organism(name=name, taxonomy_id=taxonomy_id, is_scientific_name=True)
        organism.save()
        qn_result = ComputationalResult()
        qn_result.save()
        ComputationalResultAnnotation.objects.create(
            result=qn_result,
            data={"is_qn": True, "organism_id": organism.id},
        )
        organism.qn_target = qn_result
        organism.save()
        return organism

    @classmethod
    def _make_sample(cls, accession, organism, platform_name, platform_accession):
        sample = Sample()
        sample.accession_code = accession
        sample.title = accession
        sample.is_processed = True
        sample.organism = organism
        sample.platform_name = platform_name
        sample.platform_accession_code = platform_accession
        sample.save()
        return sample

    @classmethod
    def _attach(cls, experiment, sample):
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)
        SampleResultAssociation.objects.create(sample=sample, result=experiment._parity_qn)

    @classmethod
    def setUpClass(cls):
        super(ESTestCases, cls).setUpClass()  # ref https://stackoverflow.com/a/29655301/763705

        """
        #Set up class.
        """
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
            name=ECOLI_STRAIN_NAME,
            taxonomy_id=879462,
            is_scientific_name=True,
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

        processor_job = ProcessorJob(downloader_job=downloader_job)
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

        # Extended fixture: three experiments with distinct properties so
        # compound filters, range boundaries, ordering, and facet cardinality
        # metrics are all exercised with distinguishable documents.
        human = cls._make_qn_organism("Homo sapiens", 9606)
        mouse = cls._make_qn_organism(MOUSE, 10090)

        exp_a = Experiment()
        exp_a.accession_code = cls.EXP_A
        exp_a.title = "Breast cancer gene expression"
        exp_a.description = "Routine assay"
        exp_a.technology = "MICROARRAY"
        exp_a.has_publication = True
        exp_a.publication_title = "Breast cancer findings"
        exp_a.num_processed_samples = 10
        exp_a.num_total_samples = 10
        exp_a.num_downloadable_samples = 10
        exp_a.save()
        exp_a._parity_qn = ComputationalResult.objects.create()
        for i in range(2):
            cls._attach(
                exp_a,
                cls._make_sample(f"PA1-S{i}", human, "Affymetrix HG-U133 Plus 2", "hgu133plus2"),
            )
        exp_a.update_platform_names()
        exp_a.save()

        exp_b = Experiment()
        exp_b.accession_code = cls.EXP_B
        exp_b.title = "Routine transcriptome assay"
        exp_b.description = "Cancer treatment response survey"
        exp_b.technology = "RNA-SEQ"
        exp_b.has_publication = False
        exp_b.num_processed_samples = 5
        exp_b.num_total_samples = 5
        exp_b.num_downloadable_samples = 5
        exp_b.save()
        exp_b._parity_qn = ComputationalResult.objects.create()
        for i in range(2):
            cls._attach(
                exp_b, cls._make_sample(f"PA2-S{i}", human, "Illumina HiSeq 2500", "illuminaHiSeq")
            )
        exp_b.update_platform_names()
        exp_b.save()

        exp_c = Experiment()
        exp_c.accession_code = cls.EXP_C
        exp_c.title = "Mouse brain transcriptome"
        exp_c.description = "Unrelated phrase"
        exp_c.technology = "RNA-SEQ"
        exp_c.has_publication = True
        exp_c.publication_title = "Mouse brain atlas"
        exp_c.num_processed_samples = 20
        exp_c.num_total_samples = 20
        exp_c.num_downloadable_samples = 20
        exp_c.save()
        exp_c._parity_qn = ComputationalResult.objects.create()
        for i in range(2):
            cls._attach(
                exp_c, cls._make_sample(f"PA3-S{i}", mouse, "Illumina HiSeq 2500", "illuminaHiSeq")
            )
        exp_c.update_platform_names()
        exp_c.save()

        # clear default cache and reindex
        # otherwise the organisms with qn_targes will be cached.
        cache.clear()
        call_command("opensearch", "index", "rebuild", "--force")
        call_command("opensearch", "document", "index", "--force", "--refresh")

    def test_es_endpoint(self):
        """
        #Test basic ES functionality
        """
        response = self.client.get(reverse("search", kwargs={"version": API_VERSION}))

        """
        #Test basic ES functionality
        """
        es_search_result = ExperimentDocument.search().filter("term", description="soda")
        es_search_result_qs = es_search_result.to_queryset()
        self.assertEqual(len(es_search_result_qs), 1)

        # Sanity
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 5)

        # test sample counts in filters (cardinality of downloadable_samples per bucket)
        self.assertEqual(response.json()["facets"]["has_publication"]["false"], 3)
        self.assertEqual(response.json()["facets"]["has_publication"]["true"], 4)
        self.assertEqual(response.json()["facets"]["technology"]["microarray"], 3)
        self.assertEqual(response.json()["facets"]["technology"]["rna-seq"], 4)
        self.assertEqual(
            set(response.json()["facets"]["downloadable_organism_names"].keys()),
            {ECOLI_STRAIN_NAME, "Homo sapiens", "Mus musculus"},
        )
        self.assertEqual(
            response.json()["facets"]["downloadable_organism_names"][ECOLI_STRAIN_NAME], 1
        )
        self.assertEqual(
            response.json()["facets"]["downloadable_organism_names"]["Homo sapiens"], 4
        )
        self.assertEqual(response.json()["facets"]["downloadable_organism_names"][MOUSE], 2)

        # Basic Search
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}), {"search": "soda"}
        )
        self.assertEqual(response.json()["count"], 1)
        cache.clear()

        # Positive filter result
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "soda", "technology": "microarray"},
        )
        self.assertEqual(response.json()["count"], 1)
        cache.clear()

        # Negative filter result
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "soda", "technology": "rna"},
        )
        self.assertEqual(response.json()["count"], 0)

        # Filter based on organism name
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "soda", "organism": ECOLI_STRAIN_NAME},
        )
        self.assertEqual(response.json()["count"], 1)

    def test_es_endpoint_post(self):
        # Basic filter
        search = {"accession_code": "GSE123-X"}
        response = self.client.post(
            reverse("search", kwargs={"version": API_VERSION}),
            json.dumps(search),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 1)
        cache.clear()

        # __in filter
        search = {"accession_code__in": ["GSE123-X"]}
        response = self.client.post(
            reverse("search", kwargs={"version": API_VERSION}),
            json.dumps(search),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 1)
        cache.clear()

        # Numeric filter — matches GSE123-X (1), EXP_A (10), EXP_B (5), EXP_C (20).
        search = {"num_downloadable_samples__gt": 0}
        response = self.client.post(
            reverse("search", kwargs={"version": API_VERSION}),
            json.dumps(search),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 4)
        cache.clear()

        # Large query
        search = {"accession_code__in": ["GSE123-X" for _ in range(10000)]}
        with self.assertRaises(TooManyFieldsSent):
            # This should fail for GET requests because the generated URL is too large
            response = self.client.get(reverse("search", kwargs={"version": API_VERSION}), search)

        response = self.client.post(
            reverse("search", kwargs={"version": API_VERSION}),
            json.dumps(search),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["count"], 1)

    # ---------------------------------------------------------------
    # Filter stacking
    # ---------------------------------------------------------------

    def test_compound_filter_without_search_text(self):
        """Stacked filters work with no `search=` term — exercises the GET
        FilteringFilterBackend without CompoundSearchFilterBackend interference."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"technology": "rna-seq", "has_publication": "true"},
        )
        self.assertEqual(response.status_code, 200)
        codes = {r["accession_code"] for r in response.json()["results"]}
        self.assertEqual(codes, {self.EXP_C})

    def test_compound_filter_post(self):
        """Same compound filter via POST exercises POSTFilteringFilterBackend's
        MockRequest adapter — a project-local shim the OpenSearch rewrite must
        replicate exactly."""
        response = self.client.post(
            reverse("search", kwargs={"version": API_VERSION}),
            json.dumps({"technology": "rna-seq", "organism": MOUSE}),
            content_type="application/json",
        )
        self.assertEqual(response.status_code, 200)
        codes = {r["accession_code"] for r in response.json()["results"]}
        self.assertEqual(codes, {self.EXP_C})

    # ---------------------------------------------------------------
    # Range boundaries — classic off-by-one territory
    # ---------------------------------------------------------------

    def test_range_gt_is_strict(self):
        """`__gt=10` must EXCLUDE the doc whose value is exactly 10."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"num_processed_samples__gt": 10},
        )
        codes = {r["accession_code"] for r in response.json()["results"]}
        self.assertNotIn(self.EXP_A, codes)  # exactly 10 — excluded
        self.assertIn(self.EXP_C, codes)  # 20 — included

    # ---------------------------------------------------------------
    # Analyzer / boost ranking
    # ---------------------------------------------------------------

    def test_title_boost_outranks_description_match(self):
        """`search=cancer` matches EXP_A in title (boost 10) and EXP_B in
        description (boost 2). Title match must rank first. Scores are
        version-sensitive; the relative order is not."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "cancer"},
        )
        results = response.json()["results"]
        codes = [r["accession_code"] for r in results]
        self.assertIn(self.EXP_A, codes)
        self.assertIn(self.EXP_B, codes)
        self.assertLess(codes.index(self.EXP_A), codes.index(self.EXP_B))

    # ---------------------------------------------------------------
    # Ordering and pagination
    # ---------------------------------------------------------------

    def test_explicit_ordering_overrides_default(self):
        """`ordering=-num_downloadable_samples` puts EXP_C (20) first across the
        full parity-fixture set."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"ordering": "-num_downloadable_samples"},
        )
        results = response.json()["results"]
        self.assertEqual(results[0]["accession_code"], self.EXP_C)
        self.assertEqual(results[0]["num_downloadable_samples"], 20)

    def test_pagination_limit_offset(self):
        """`limit=1&offset=1` returns the second item of the ordered page;
        `count` still reflects the total (5)."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"ordering": "-num_downloadable_samples", "limit": 1, "offset": 1},
        )
        body = response.json()
        self.assertEqual(body["count"], 5)
        self.assertEqual(len(body["results"]), 1)
        # Ranked by num_downloadable_samples desc: EXP_C (20), EXP_A (10), ...
        self.assertEqual(body["results"][0]["accession_code"], self.EXP_A)

    # ---------------------------------------------------------------
    # Facets — the custom FacetedSearchFilterBackendExtended cardinality metric
    # ---------------------------------------------------------------

    def test_facets_present_on_empty_result(self):
        """A zero-hit query must still produce the facets envelope with the
        expected keys. OpenSearch's agg response shape matches ES but verify
        explicitly so a client rewrite cannot quietly drop the key."""
        response = self.client.get(
            reverse("search", kwargs={"version": API_VERSION}),
            {"search": "xyzzynonexistentterm"},
        )
        body = response.json()
        self.assertEqual(body["count"], 0)
        for key in ("technology", "has_publication", "downloadable_organism_names"):
            self.assertIn(key, body["facets"])

    def test_facet_cardinality_metric_is_integer(self):
        """`FacetedSearchFilterBackendExtended` adds a `cardinality` sub-metric
        (`total_samples`) to every facet bucket. After `transform_es_facets`,
        each bucket value should be a non-negative integer regardless of
        engine version."""
        response = self.client.get(reverse("search", kwargs={"version": API_VERSION}))
        for field in ("technology", "has_publication", "downloadable_organism_names"):
            for bucket_value in response.json()["facets"][field].values():
                self.assertIsInstance(bucket_value, int)
                self.assertGreaterEqual(bucket_value, 0)
