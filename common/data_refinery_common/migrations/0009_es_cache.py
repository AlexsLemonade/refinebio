# Generated by Django 2.1.2 on 2019-01-08 20:20
# Modified by Rich

import django.contrib.postgres.fields
from django.db import migrations, models


def update_cached_values(apps, schema_editor):
    """ """

    Experiment = apps.get_model("data_refinery_common", "Experiment")
    for experiment in Experiment.objects.all():
        # Model methods can't be used during migrations. :(
        # https://stackoverflow.com/a/37685925/1135467

        # via experiment.update_num_samples
        experiment.num_total_samples = experiment.samples.count()
        experiment.num_processed_samples = experiment.samples.filter(is_processed=True).count()

        # via experiment.update_sample_metadata_fields
        fields = []
        possible_fields = [
            "sex",
            "age",
            "specimen_part",
            "genotype",
            "disease",
            "disease_stage",
            "cell_line",
            "treatment",
            "race",
            "subject",
            "compound",
            "time",
        ]
        samples = experiment.samples.all()
        for field in possible_fields:
            for sample in samples:
                if getattr(sample, field) != None and getattr(sample, field) != "":
                    fields.append(field)
                    break

        experiment.sample_metadata_fields = fields

        # via experiment.update_organism_names
        experiment.organism_names = list(
            set([organism.name for organism in experiment.organisms.all()])
        )

        # experiment.update_platform_names
        experiment.platform_names = list(
            set([sample.platform_name for sample in experiment.samples.all()])
        )
        experiment.save()


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0008_auto_20190104_2135"),
    ]

    operations = [
        migrations.AddField(
            model_name="experiment",
            name="num_processed_samples",
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name="experiment",
            name="num_total_samples",
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name="experiment",
            name="organism_names",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.TextField(), default=list, size=None
            ),
        ),
        migrations.AddField(
            model_name="experiment",
            name="sample_metadata_fields",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.TextField(), default=list, size=None
            ),
        ),
        migrations.AddField(
            model_name="experiment",
            name="platform_names",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.TextField(), default=list, size=None
            ),
        ),
        migrations.AddField(
            model_name="originalfile",
            name="downloader_jobs",
            field=models.ManyToManyField(
                through="data_refinery_common.DownloaderJobOriginalFileAssociation",
                to="data_refinery_common.DownloaderJob",
            ),
        ),
        migrations.AddField(
            model_name="originalfile",
            name="processor_jobs",
            field=models.ManyToManyField(
                through="data_refinery_common.ProcessorJobOriginalFileAssociation",
                to="data_refinery_common.ProcessorJob",
            ),
        ),
        migrations.AddField(
            model_name="sample",
            name="experiments",
            field=models.ManyToManyField(
                through="data_refinery_common.ExperimentSampleAssociation",
                to="data_refinery_common.Experiment",
            ),
        ),
        migrations.RunPython(update_cached_values),
    ]
