from django.test import TestCase, tag


class ImportTestCase(TestCase):
    @tag("downloaders")
    def test_downloader_imports(self):
        # Make sure we can import the downloader tests.
        from tests.downloaders import (
            test_array_express,
            test_geo,
            test_sra,
            test_transcriptome_index,
        )
