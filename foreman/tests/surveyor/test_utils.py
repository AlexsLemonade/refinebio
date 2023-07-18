from django.test import TestCase

from data_refinery_foreman.surveyor import utils


class UtilsTestCase(TestCase):
    def test_flatten_dict(self):
        """
        Tests that flattening dictionaries and lists behaves as expected.
        """
        source_dict = {
            "first foo": "first_value",
            "second bar": {
                "nested": True,
            },
            "third fuzz": ["another_nested", {"nested_object": True}],
        }

        flattened_dict = utils.flatten_dict(source_dict)

        self.assertEqual(flattened_dict["first_foo"], "first_value")
        self.assertEqual(flattened_dict["second_bar_nested"], True)
        self.assertEqual(flattened_dict["third_fuzz_0"], "another_nested")
        self.assertEqual(flattened_dict["third_fuzz_1_nested_object"], True)

        prefix_flattened_dict = utils.flatten_dict(source_dict, "test")

        self.assertEqual(prefix_flattened_dict["test_first_foo"], "first_value")
        self.assertEqual(prefix_flattened_dict["test_second_bar_nested"], True)
        self.assertEqual(prefix_flattened_dict["test_third_fuzz_0"], "another_nested")
        self.assertEqual(prefix_flattened_dict["test_third_fuzz_1_nested_object"], True)
