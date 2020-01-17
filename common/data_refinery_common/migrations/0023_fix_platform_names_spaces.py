from django.db import migrations

from data_refinery_common.utils import get_supported_rnaseq_platforms


def make_sample_platform_names_readable(apps, schema_editor):
    Sample = apps.get_model("data_refinery_common", "Sample")

    for platform_name in get_supported_rnaseq_platforms():
        munged_platform = platform_name.replace(" ", "")
        Sample.objects.all().filter(platform_name=munged_platform).update(
            platform_name=platform_name
        )
        print("Updating platform name from '%s' to '%s'" % (munged_platform, platform_name))


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0022_auto_20190607_1505"),
    ]

    operations = [migrations.RunPython(make_sample_platform_names_readable)]
