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

    def test_get_first_dict(self):
        """
        Tests that you can fetch the first dict from an iterable based on a key.
        """

        iterable_of_dicts = [
            {"test_key": 1},
            {"test_key": 2},
            {"test_key": 3},
            {"test_key": 4},
        ]

        # Test for successful match.
        matched_dict = utils.find_first_dict("test_key", 1, iterable_of_dicts)
        self.assertEqual(matched_dict["test_key"], 1)

        # Test when key exists, but no match.
        good_key_no_match = utils.find_first_dict("test_key", 5, iterable_of_dicts)
        self.assertEqual(good_key_no_match, None)

        # Test when key does not existing.
        no_dict = utils.find_first_dict("another_key", 1, iterable_of_dicts)
        self.assertEqual(no_dict, None)
