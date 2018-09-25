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
