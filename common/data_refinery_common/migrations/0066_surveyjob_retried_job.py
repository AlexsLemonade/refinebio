# Generated by Django 2.2.13 on 2021-04-01 15:04

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0065_auto_20210318_1528"),
    ]

    operations = [
        migrations.AddField(
            model_name="surveyjob",
            name="retried_job",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to="data_refinery_common.SurveyJob",
            ),
        ),
    ]