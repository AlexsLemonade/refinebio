from django.test import TestCase, tag


class ImportTestCase(TestCase):
    def test_imports(self):
        # Make sure we can import the common tests
        import data_refinery_common.models.test_jobs
        import data_refinery_common.models.test_models
        import data_refinery_common.models.test_ontology_term
        import data_refinery_common.models.test_organisms
        import data_refinery_common.test_job_management
        import data_refinery_common.test_microarray
        import data_refinery_common.test_utils
