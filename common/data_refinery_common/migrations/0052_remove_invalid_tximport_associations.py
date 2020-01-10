from django.db import migrations, models
from django.db.models import Count, F


def remove_tximport_invalid_associations(apps, schema_editor):
    """ We were associating tximport with all samples in an experiment
    even though some of them were unprocessed.
    Ref https://github.com/AlexsLemonade/refinebio/issues/2054 """
    SampleResultAssociation = apps.get_model("data_refinery_common", "SampleResultAssociation")

    SampleResultAssociation.objects.filter(
        sample__is_processed=False, result__processor__name__iexact="tximport"
    ).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0051_remove_dataset_email_sent"),
    ]

    operations = [
        migrations.RunPython(remove_tximport_invalid_associations),
    ]
