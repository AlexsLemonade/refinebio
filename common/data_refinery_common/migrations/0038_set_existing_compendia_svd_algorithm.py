from django.db import migrations, models
from django.core.paginator import Paginator

def set_existing_compendia_svd_algorithm(apps, schema_editor):
    """ Fixes existing compendia svd algorithm - defaulted NONE corrects to ARPACK """ 
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')

    paginator = Paginator(ComputedFile.objects.all().order_by('id').filter(is_compendia=True), 1000)
    for page_idx in range(1, paginator.num_pages):
        for computedFile in paginator.page(page_idx).object_list:
            computedFile.update(svd_algorithm="ARPACK")


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0037_auto_20190924_2305'),
    ]

    operations = [
        migrations.RunPython(set_existing_compendia_svd_algorithm),
    ]
