from django.core.management.base import BaseCommand
from django.db import transaction

# from data_refinery_common.models.volume import Volume


class Command(BaseCommand):
    def handle(self, *args, **options):
        with transaction.atomic():
            volume = Volume.inactive_objects.order_by("-last_modified").first()
            volume_id = volume.id
            volume.is_active = True
            volume.save()

        print("INSTANCE_ID={}".format(volume_id))
