# Generated by Django 2.2.13 on 2021-03-18 15:28

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0064_auto_20210225_1622"),
    ]

    operations = [
        migrations.RenameField(
            model_name="downloaderjob", old_name="nomad_job_id", new_name="batch_job_id",
        ),
        migrations.RenameField(
            model_name="processorjob", old_name="nomad_job_id", new_name="batch_job_id",
        ),
        migrations.RenameField(
            model_name="surveyjob", old_name="nomad_job_id", new_name="batch_job_id",
        ),
    ]
