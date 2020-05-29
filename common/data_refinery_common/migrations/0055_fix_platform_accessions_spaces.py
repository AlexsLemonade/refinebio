from django.db import migrations

from data_refinery_common.utils import get_supported_rnaseq_platforms


def remove_spaces_from_platform_accessions(apps, schema_editor):
    Sample = apps.get_model("data_refinery_common", "Sample")

    for bad_accession in get_supported_rnaseq_platforms():
        platform_accession = bad_accession.replace(" ", "")
        bad_samples = Sample.objects.all().filter(platform_accession_code=bad_accession)

        if not bad_samples:
            continue

        bad_samples.update(platform_accession_code=platform_accession)
        print("Updating platform accession from '%s' to '%s'" % (bad_accession, platform_accession))


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0054_datasetannotation"),
    ]

    operations = [migrations.RunPython(remove_spaces_from_platform_accessions)]
