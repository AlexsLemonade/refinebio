from django.db import migrations, models

def set_existing_compendia_svd_algorithm(apps, schema_editor):
    """ Fixes affymetrix samples that have their manufacturer set to "AFFYMETRTIX" or "NEXTSEQ" """
    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    ComputedFile.objects.all().filter(is_compendia=True).update(svd_algorithm="arpack")


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0036_auto_20190919_2306'),
    ]

    operations = [
        migrations.RunPython(set_existing_compendia_svd_algorithm),
    ]
