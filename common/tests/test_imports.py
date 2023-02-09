from django.test import TestCase


class ImportTestCase(TestCase):
    def test_imports(self):
        # Make sure we can import the common tests
        import tests.models.test_jobs
        import tests.models.test_models
        import tests.models.test_ontology_term
        import tests.models.test_organisms
        import tests.test_job_management
        import tests.test_microarray
        import tests.test_utils
