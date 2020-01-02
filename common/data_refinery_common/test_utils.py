from unittest.mock import Mock, patch
from django.test import TestCase
from data_refinery_common import utils
from data_refinery_common.models import Pipeline


class UtilsTestCase(TestCase):
    @patch("data_refinery_common.utils.requests.get")
    def test_get_instance_id_cloud(self, mock_get):
        """Test that a request is made and the global value is stored"""
        # Ensure utils.INSTANCE_ID hasn't been set yet in case the
        # order the tests are run in ever changes
        utils.INSTANCE_ID = None
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = "instance_id"

        with self.settings(RUNNING_IN_CLOUD=True):
            self.assertEqual(utils.get_instance_id(), "instance_id")

        # Ensure that the second call uses the now-set global value.
        # (By resetting the mocks, calling it again, and checking that
        # the values didn't need to be set again).
        mock_get.reset_mock()
        utils.get_instance_id()
        mock_get.assert_not_called()

    def test_get_instance_id_local(self):
        """Test that local is used for instance id."""
        # Ensure utils.INSTANCE_ID hasn't been set yet in case the
        # order the tests are run in ever changes
        utils.INSTANCE_ID = None

        with self.settings(RUNNING_IN_CLOUD=False):
            self.assertEqual(utils.get_instance_id(), "local")

        # Ensure that the second call uses the now-set global value
        # by changing what settings would tell it.
        with self.settings(RUNNING_IN_CLOUD=True):
            self.assertEqual(utils.get_instance_id(), "local")

    def test_supported_microarray_platforms(self):
        """Test that supported microarray platforms setting is set correctly."""
        supported_microarray_platforms = utils.get_supported_microarray_platforms()

        has_equgene11st = False
        has_A_AFFY_59 = False
        has_GPL23026 = False
        has_AGEOD23026 = False
        for platform in supported_microarray_platforms:
            if platform["platform_accession"] == "equgene11st" and platform["is_brainarray"]:
                has_equgene11st = True

            if platform["external_accession"] == "A-AFFY-59" and not platform["is_brainarray"]:
                has_A_AFFY_59 = True

            if platform["external_accession"] == "GPL23026" and not platform["is_brainarray"]:
                has_GPL23026 = True

            if platform["external_accession"] == "A-GEOD-23026" and not platform["is_brainarray"]:
                has_AGEOD23026 = True

        self.assertTrue(has_equgene11st)
        self.assertTrue(has_A_AFFY_59)
        self.assertTrue(has_GPL23026)
        self.assertTrue(has_AGEOD23026)

    def test_get_internal_microarray_accession(self):
        """Test that supported microarray platforms setting is set correctly."""

        self.assertEqual(utils.get_internal_microarray_accession("hgu133a"), "hgu133a")
        self.assertEqual(utils.get_internal_microarray_accession("A-AFFY-59"), "soybean")
        self.assertEqual(
            utils.get_internal_microarray_accession("GPL23026"), "Illumina_HumanHT-12_V4.0"
        )

    def test_supported_rnaseq_platforms(self):
        """Test that supported RNASeq platforms setting is set correctly."""
        self.assertTrue("Illumina HiSeq 1000" in utils.get_supported_rnaseq_platforms())

    def test_readable_affymetrix_names(self):
        """Test that the setting for Affymetrix accessions to
        human readable names is set correctly."""
        readable_platform_names = utils.get_readable_affymetrix_names()
        expected_readable_name = "[ChiGene-1_0-st] Affymetrix Chicken Gene 1.0 ST Array"
        self.assertTrue(readable_platform_names["chigene10st"] == expected_readable_name)
        expected_readable_name = "[Xenopus_laevis] Affymetrix Xenopus laevis Genome Array"
        self.assertTrue(readable_platform_names["xenopuslaevis"] == expected_readable_name)

    def test_get_normalized_platform(self):
        """ Test a particular normaization we need to perform """

        self.assertEqual(utils.get_normalized_platform("hugene10stv1"), "hugene10st")
        self.assertEqual(utils.get_normalized_platform("hugene10stv2"), "hugene10st")
        self.assertEqual(utils.get_normalized_platform("stv1hugene10"), "stv1hugene10")

    def test_volume_index(self):
        """Test that supported RNASeq platforms setting is set correctly."""
        self.assertEqual(utils.get_volume_index(), "0")
        with open("/tmp/VOLUME_INDEX", "wb") as f:
            f.write("123".encode())
        self.assertEqual(utils.get_volume_index(path="/tmp/VOLUME_INDEX"), "123")

    def test_load_blacklist(self):
        blacklist = utils.load_blacklist()
        self.assertEqual(len(blacklist), 239449)

    def test_queryset_iterator(self):
        """Test that the queryset iterator by using it to actually iterate over a queryset.

        Uses Pipeline just because it's easy to init."""
        # Page size defaults to 2000, so use something bigger than
        # that so there's more than one page.
        for i in range(3000):
            Pipeline(name=str(i)).save()

        pipelines = Pipeline.objects.all()

        # Build a list of the names just to do something with the data
        # so we know the query actually resolved.
        names = []
        for pipeline in utils.queryset_iterator(pipelines):
            names.append(pipeline.name)

        self.assertEqual(len(names), 3000)
