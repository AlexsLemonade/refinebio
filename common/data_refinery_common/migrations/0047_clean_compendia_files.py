from django.db import migrations
from django.db.models import Count


def clean_compendium_files(apps, schema_editor):
    """ We used to have more than one computed file for a compendium
    Now we just want to keep the zipped file that lives on s3."""
    CompendiumResult = apps.get_model("data_refinery_common", "CompendiumResult")
    ComputedFile = apps.get_model("data_refinery_common", "ComputedFile")

    # normalized compendia
    compendium_results = (
        CompendiumResult.objects.annotate(file_count=Count("result__computedfile"))
        .filter(quant_sf_only=False, file_count__gt=1)
        .values("result_id")
    )

    ComputedFile.objects.filter(result_id__in=compendium_results).filter(
        is_compendia=False
    ).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0046_compendium_version_number"),
    ]

    operations = [
        migrations.RunPython(clean_compendium_files),
    ]
