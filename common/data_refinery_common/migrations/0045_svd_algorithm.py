from django.db import migrations, models


def svd_algorithm_fix(apps, schema_editor):
    """ Fetch existing compendia and create compendium results"""
    apps.get_model("data_refinery_common", "CompendiumResult")
    ComputedFile = apps.get_model("data_refinery_common", "ComputedFile")

    compendium_computed_files = ComputedFile.objects.filter(
        is_compendia=True, result__compendium_result__isnull=False
    ).prefetch_related("result__compendium_result")

    for computed_file in compendium_computed_files:
        if computed_file.svd_algorithm is not "None":
            compendium_result = computed_file.result.compendium_result.first()
            compendium_result.svd_algorithm = computed_file.svd_algorithm
            compendium_result.save()


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0044_compendium_results"),
    ]

    operations = [
        migrations.RunPython(svd_algorithm_fix),
    ]
