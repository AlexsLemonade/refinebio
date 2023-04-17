from django.test import TestCase

from data_refinery_common import microarray


class MicroarrayTestCase(TestCase):
    def test_get_platform_from_CEL(self):
        self.assertEqual(
            "hgu95av2", microarray.get_platform_from_CEL("tests/data/microarray/C30057.CEL.gz")
        )
        self.assertEqual(
            "rgu34a", microarray.get_platform_from_CEL("tests/data/microarray/SG2_u34a.CEL.gz")
        )
        self.assertEqual(
            "mouse4302",
            microarray.get_platform_from_CEL("tests/data/microarray/97_(Mouse430_2).CEL.gz"),
        )
        self.assertEqual(
            "zebgene11st", microarray.get_platform_from_CEL("tests/data/microarray/CONTROL6.cel.gz")
        )
