from django.db import migrations, models


def batch_update(queryset, batch=1000, **changes):
    """ Update per batch """
    is_done = False
    current = 0

    while not is_done:
        current = current + 1
        start = (current - 1) * batch
        end = start + batch
        pks = queryset.values_list("pk", flat=True)[start:end][::1]
        if pks:
            queryset.model.objects.filter(pk__in=pks).update(**changes)

        is_done = len(pks) < batch


def set_computed_file_svd_algorithm_arpack(apps, schema_editor):
    """ Set svd_algorithm for computed files - ARPACK when is_compendia"""
    ComputedFile = apps.get_model("data_refinery_common", "ComputedFile")
    arpack_queryset = ComputedFile.objects.filter(is_compendia=True).order_by("id")
    batch_update(arpack_queryset, svd_algorithm="ARPACK")


def set_computed_file_svd_algorithm_none(apps, schema_editor):
    """ Set svd_algorithm for computed files - None when not is_compendia"""
    ComputedFile = apps.get_model("data_refinery_common", "ComputedFile")
    none_queryset = ComputedFile.objects.filter(is_compendia=False).order_by("id")
    batch_update(none_queryset, svd_algorithm="NONE")


def set_dataset_svd_algorithm_none(apps, schema_editor):
    """ Apply default svd_algorithm via pagination to dataset """
    Dataset = apps.get_model("data_refinery_common", "Dataset")
    queryset = Dataset.objects.all().order_by("id")
    batch_update(queryset, svd_algorithm="NONE")


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0037_auto_20190925_1309"),
    ]

    operations = [
        migrations.RunPython(set_computed_file_svd_algorithm_arpack),
        migrations.RunPython(set_computed_file_svd_algorithm_none),
        migrations.RunPython(set_dataset_svd_algorithm_none),
    ]
