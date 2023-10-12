from django.test import TestCase

from data_refinery_common.models import Sample


class TestSample(TestCase):
    def setUp(self):
        self.sample = Sample()
        self.sample.save()

    def test_to_metadata_dict(self):
        self.assertListEqual(
            sorted(self.sample.to_metadata_dict().keys()),
            [
                "refinebio_accession_code",
                "refinebio_age",
                "refinebio_annotations",
                "refinebio_cell_line",
                "refinebio_compound",
                "refinebio_disease",
                "refinebio_disease_stage",
                "refinebio_genetic_information",
                "refinebio_organism",
                "refinebio_platform",
                "refinebio_processed",
                "refinebio_race",
                "refinebio_sex",
                "refinebio_source_archive_url",
                "refinebio_source_database",
                "refinebio_specimen_part",
                "refinebio_subject",
                "refinebio_time",
                "refinebio_title",
                "refinebio_treatment",
            ],
        )
