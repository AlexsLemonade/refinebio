from django.db import migrations, models
from django.core.paginator import Paginator

def set_computed_file_svd_algorithm_none(apps, schema_editor):
    """ Set svd_algorithm to NONE for existing computed files that are not compendia"""    
    ### is_compendia=False
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    paginator = Paginator(ComputedFile.objects.all().order_by('id').filter(is_compendia=False), 1000)
    for page_idx in range(1, paginator.num_pages):
        for computedFile in paginator.page(page_idx).object_list:
            computedFile.update(svd_algorithm="NONE")

def set_computed_file_svd_algorithm_arpack(apps, schema_editor):
    """ Set svd_algorithm to ARPACK for exiting compendia """
    ### is_compendia=True   
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    paginator = Paginator(ComputedFile.objects.all().order_by('id').filter(is_compendia=True), 1000)
    for page_idx in range(1, paginator.num_pages):
        for computedFile in paginator.page(page_idx).object_list:
            computedFile.update(svd_algorithm="ARPACK")

def set_dataset_svd_algorithm_none(apps, schema_editor):
    """ Apply default svd_algorithm via pagination to dataset """    
    Dataset = apps.get_model('data_refinery_common', 'Dataset')
    paginator = Paginator(Dataset.objects.all().order_by('id'), 1000)
    for page_idx in range(1, paginator.num_pages):
        for dataset in paginator.page(page_idx).object_list:
            dataset.update(svd_algorithm="NONE")

class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0037_auto_20190925_1309'),
    ]

    operations = [
        migrations.RunPython(set_computed_file_svd_algorithm_none),
        migrations.RunPython(set_computed_file_svd_algorithm_arpack),
        migrations.RunPython(set_dataset_svd_algorithm_none),
    ]
