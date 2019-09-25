from django.db import migrations, models
from django.core.paginator import Paginator

def set_computed_file_svd_algorithm_arpack(apps, schema_editor):
    """ Set svd_algorithm for computed files - ARPACK when is_compendia"""
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    paginator = Paginator(ComputedFile.objects.filter(is_compendia=True).order_by('id'), 1000)
    for page_idx in range(1, paginator.num_pages):
        paginator.page(page_idx).object_list.update(svd_algorithm="ARPACK")

def set_computed_file_svd_algorithm_none(apps, schema_editor):
    """ Set svd_algorithm for computed files - None when not is_compendia"""
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    paginator = Paginator(ComputedFile.objects.filter(is_compendia=False).order_by('id'), 1000)
    for page_idx in range(1, paginator.num_pages):
        paginator.page(page_idx).object_list.update(svd_algorithm="None")

def set_dataset_svd_algorithm_none(apps, schema_editor):
    """ Apply default svd_algorithm via pagination to dataset """
    Dataset = apps.get_model('data_refinery_common', 'Dataset')
    paginator = Paginator(Dataset.objects.all().order_by('id'), 1000)
    for page_idx in range(1, paginator.num_pages):
        paginator.page(page_idx).object_list.update(svd_algorithm="NONE")

class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0037_auto_20190925_1309'),
    ]

    operations = [
        migrations.RunPython(set_computed_file_svd_algorithm_arpack),
        migrations.RunPython(set_computed_file_svd_algorithm_none),
        migrations.RunPython(set_dataset_svd_algorithm_none),
    ]
