from django.test import TestCase


class ImportTestCase(TestCase):
    def test_imports(self):
        # Make sure we can import the api tests
        import tests.management.commands.test_post_downloads_summary
        import tests.views.test_api_general
        import tests.views.test_compendia
        import tests.views.test_dataset
        import tests.views.test_dataset_stats
        import tests.views.test_processor
        import tests.views.test_qn_target
        import tests.views.test_search
        import tests.test_sentry_middleware
        import tests.views.test_stats
