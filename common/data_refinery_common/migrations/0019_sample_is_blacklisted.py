# Generated by Django 2.1.5 on 2019-04-08 14:34

from django.core.paginator import Paginator
from django.db import migrations, models

from data_refinery_common.utils import load_blacklist


def update_blacklisted_samples(apps, schema_editor):
    """
    Uses the SRA-supplied blacklist to mark samples as is_blacklisted or not.
    """
    Sample = apps.get_model("data_refinery_common", "Sample")

    blacklist = load_blacklist()

    paginator = Paginator(Sample.objects.all().order_by("id"), 1000)
    for page_idx in range(1, paginator.num_pages):
        blacklisted = []
        for sample in paginator.page(page_idx).object_list:
            if sample.accession_code in blacklist:
                sample.is_blacklisted = True
                sample.save()
                blacklisted.append(str(sample.accession_code))
        print("Blacklisted page " + str(page_idx) + ": " + str(blacklisted))


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0018_downloaderjob_ram_amount"),
    ]

    operations = [
        migrations.AddField(
            model_name="sample",
            name="is_blacklisted",
            field=models.BooleanField(default=False),
        ),
        migrations.RunPython(update_blacklisted_samples),
    ]
