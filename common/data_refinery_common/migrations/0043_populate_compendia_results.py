from django.db import migrations, models

def populate_compendia_results(apps, schema):
    """ Populate compendia_results with computed files """

    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    CompendiaResult = apps.get_model('data_refinery_common', 'CompendiaResult')

    # get all existing compendia computed files
    compendia_computed_files = ComputedFile.objects.filter(is_compendia=True)

    # create as compendia result from each of the results
    compendia_results = []
    for computed_file in compendia_computed_files:
        compendia_result = CompendiaResult(quant_sf_only=computed_file.quant_sf_only,
                                           compendia_version=computed_file.compendia_version,
                                           svd_algorithm=computed_file.svd_algorithm)

        # relationships
        compendia_result.result = computed_file.result
        compendia_result.primary_organism = computed_file.compendia_organism
        compendia_result.organisms.add(computed_file.compendia_organism)
        compendia_results.append(compendia_result)

    CompendiaResult.objects.bulk_create(compendia_results)

class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0042_auto_20191031_1953'),
    ]

    operations = [
        migrations.RunPython(populate_compendia_results)
    ]
