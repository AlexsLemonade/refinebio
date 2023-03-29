from django.test import TestCase, tag


class ImportTestCase(TestCase):
    @tag("downloaders")
    def test_downloader_imports(self):
        # Make sure we can import the downloader tests
        import data_refinery_workers.downloaders.test_array_express
        import data_refinery_workers.downloaders.test_geo
        import data_refinery_workers.downloaders.test_sra
        import data_refinery_workers.downloaders.test_transcriptome_index
