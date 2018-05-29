from unittest.mock import Mock, patch
from django.test import TestCase
from data_refinery_common import microarray

CEL_FILE_HUMAN = "test-files/C30057.CEL"
CEL_FILE_RAT = "test-files/SG2_u34a.CEL"
CEL_FILE_MOUSE = "test-files/97_(Mouse430_2).CEL"
CEL_FILE_ZEBRAFISH = "test-files/CONTROL6.cel"


class MicroarrayTestCase(TestCase):

    def test_get_platform_from_CEL(self):
        self.assertEqual("hgu95ahsensgprobe", microarray.get_platform_from_CEL(CEL_FILE_HUMAN))
        self.assertEqual("rgu34arnensgprobe", microarray.get_platform_from_CEL(CEL_FILE_RAT))
        self.assertEqual("", microarray.get_platform_from_CEL(CEL_FILE_MOUSE))
        self.assertEqual("rgu34arnensgprobe", microarray.get_platform_from_CEL(CEL_FILE_ZEBRAFISH))
