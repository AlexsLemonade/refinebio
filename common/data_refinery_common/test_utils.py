from unittest.mock import Mock, patch
from django.test import TestCase
from data_refinery_common import utils


class UtilsTestCase(TestCase):

    @patch('data_refinery_common.utils.get_env_variable')
    @patch('data_refinery_common.utils.requests.get')
    def test_get_worker_id_cloud(self, mock_get, mock_get_env_variable):
        """Test that a request is made and the global value is stored"""
        # Ensure utils.INSTANCE_ID hasn't been set yet in case the
        # order the tests are run in ever changes
        utils.INSTANCE_ID = None
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = "instance_id"
        mock_get_env_variable.return_value = "True"

        self.assertEqual(utils.get_worker_id(), "instance_id/MainProcess")

        # Ensure that the second call uses the now-set global value.
        # (By calling it again and checking that two function calls
        # resulted in one call to each of the invoked functions.)
        utils.get_worker_id()
        mock_get.assert_called_once()
        mock_get_env_variable.assert_called_once()

    @patch('data_refinery_common.utils.get_env_variable')
    def test_get_worker_id_local(self, mock_get_env_variable):
        """Test that local is used for instance id."""
        # Ensure utils.INSTANCE_ID hasn't been set yet in case the
        # order the tests are run in ever changes
        utils.INSTANCE_ID = None
        mock_get_env_variable.return_value = "False"

        self.assertEqual(utils.get_worker_id(), "local/MainProcess")

        # Ensure that the second call uses the now-set global value.
        # (By calling it again and checking that two function calls
        # resulted in one call to get_env_variable)
        utils.get_worker_id()
        mock_get_env_variable.assert_called_once()

    def test_supported_microarray_platforms(self):
        """Test that supported microarray platforms setting is set correctly."""
        supported_microarray_platforms = utils.get_supported_microarray_platforms()

        has_equgene11st = False
        has_A_AFFY_59 = False
        has_GPL23026 = False
        for platform in supported_microarray_platforms:
            if platform["platform_accession"] == "equgene11st" and platform["is_brainarray"]:
                has_equgene11st = True

            if platform["external_accession"] == "A-AFFY-59" and not platform["is_brainarray"]:
                has_A_AFFY_59 = True

            if platform["external_accession"] == "GPL23026" and not platform["is_brainarray"]:
                has_GPL23026 = True

        self.assertTrue(has_equgene11st)
        self.assertTrue(has_A_AFFY_59)
        self.assertTrue(has_GPL23026)

    def test_supported_rnaseq_platforms(self):
        """Test that supported RNASeq platforms setting is set correctly."""
        self.assertTrue("Illumina HiSeq 1000" in utils.get_supported_rnaseq_platforms())

    def test_readable_platform_names(self):
        """Test that the setting for mapping platform accessions to
        human readable names is set correctly."""
        readable_platform_names = utils.get_readable_platform_names()
        expected_readable_name = "[ChiGene-1_0-st] Affymetrix Chicken Gene 1.0 ST Array"
        self.assertTrue(readable_platform_names["chigene10st"] == expected_readable_name)
        expected_readable_name = "[Xenopus_laevis] Affymetrix Xenopus laevis Genome Array"
        self.assertTrue(readable_platform_names["xenopuslaevis"] == expected_readable_name)
