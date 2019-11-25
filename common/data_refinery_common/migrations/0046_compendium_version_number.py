from django.db import migrations, models
from django.db.models import F, Window
from django.db.models.functions import RowNumber

def compendium_version_number_fix(apps, schema_editor):
    """ Fetch existing compendia and fix their version_number."""
    CompendiumResult = apps.get_model('data_refinery_common', 'CompendiumResult')

    # normalized compendia
    # get each compendium against primary key sorted by created date
    compendium_results = CompendiumResult.objects.filter(quant_sf_only=False)\
                                                 .annotate(version_number=Window(
                                                     expression=RowNumber(),
                                                     partition_by=[F('primary_organism')],
                                                     order_by=[F('result__created_at')]
                                                 ))

    for compendium_result in compendium_results:
        # assign the correct version if greater than 1
        if compendium_result.version_number is not 1:
            compendium_result.compendium_version = compendium_result.version_number
            compedium_result.save()


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0045_svd_algorithm'),
    ]

    operations = [
        migrations.RunPython(compendium_version_number_fix),
    ]
