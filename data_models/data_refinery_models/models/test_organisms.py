from django.test import TestCase
from django.utils import timezone
from data_refinery_models.models.organism import *


class TimeTrackedModelTestCase(TestCase):
    def setUp(self):
        Organism.objects.create(name="HOMO SAPIENS",
                                taxonomy_id=9606,
                                is_scientific_name=True)

    def tearDown(self):
        Organism.objects.all().delete()

    def test_time_tracking_works(self):
        """created_at and updated_at are initialized upon creation"""
        organism = Organism.objects.get(name="HOMO SAPIENS")
        timedelta = timezone.now() - organism.created_at
        self.assertLess(timedelta.total_seconds(), 1)
        self.assertEqual(organism.created_at, organism.updated_at)

        # When the organism is updated and saved, updated_at changes
        organism.success = False
        organism.save()
        self.assertNotEqual(organism.created_at, organism.updated_at)
