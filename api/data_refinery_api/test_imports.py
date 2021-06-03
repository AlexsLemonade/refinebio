from django.test import TestCase, tag


class ImportTestCase(TestCase):
    def test_imports(self):
        # Make sure we can import the api tests
        import data_refinery_api.management.commands.test_post_downloads_summary
        import data_refinery_api.test.test_api_general
        import data_refinery_api.test.test_compendia
        import data_refinery_api.test.test_dataset
        import data_refinery_api.test.test_dataset_stats
        import data_refinery_api.test.test_processor
        import data_refinery_api.test.test_qn_target
        import data_refinery_api.test.test_search
        import data_refinery_api.test.test_sentry_middleware
        import data_refinery_api.test.test_stats
