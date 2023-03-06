from django.test import TestCase, tag


class ImportTestCase(TestCase):
    @tag("downloaders")
    def test_downloader_imports(self):
        # Make sure we can import the downloader tests.
        from tests.downloaders import test_array_express
        from tests.downloaders import test_geo
        from tests.downloaders import test_sra
        from tests.downloaders import test_transcriptome_index
