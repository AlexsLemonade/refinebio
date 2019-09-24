from django.db import migrations, models

def set_existing_compendia_svd_algorithm(apps, schema_editor):
    """ Fixes existing compendia svd algorithm - defaulted NONE corrects to ARPACK """ 
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    ComputedFile.objects.all().filter(is_compendia=True).update(svd_algorithm="ARPACK")


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0037_auto_20190924_2305'),
    ]

    operations = [
        migrations.RunPython(set_existing_compendia_svd_algorithm),
    ]
