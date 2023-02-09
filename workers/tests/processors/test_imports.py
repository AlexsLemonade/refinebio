from django.test import TestCase, tag


class ImportTestCase(TestCase):
    @tag("salmon")
    def test_salmon_imports(self):
        # Make sure we can import the salmon tests
        import tests.processors.test_salmon

    @tag("transcriptome")
    def test_transcriptome_imports(self):
        # Make sure we can import the transcriptome tests
        import tests.processors.test_transcriptome_index

    @tag("no_op")
    def test_no_op_imports(self):
        # Make sure we can import the no_op tests
        import tests.processors.test_no_op

    @tag("smasher")
    def test_smasher_imports(self):
        # Make sure we can import the smasher tests
        import tests.processors.test_smasher

    @tag("illumina")
    def test_illumina_imports(self):
        # Make sure we can import the illumina tests
        import tests.processors.test_illumina

    @tag("agilent")
    def test_agilent_imports(self):
        # Make sure we can import the agilent tests
        import tests.processors.test_agilent

    @tag("affymetrix")
    def test_affymetrix_imports(self):
        # Make sure we can import the affy tests
        import tests.processors.test_array_express

    @tag("qn")
    def test_qn_imports(self):
        # Make sure we can import the qn tests
        import tests.processors.test_qn_reference

    @tag("janitor")
    def test_janitor_imports(self):
        # Make sure we can import the janitor tests
        import tests.processors.test_janitor

    @tag("compendia")
    def test_compendia_imports(self):
        # Make sure we can import the compendia tests
        import tests.processors.test_compendia
        import tests.processors.test_create_quantpendia
