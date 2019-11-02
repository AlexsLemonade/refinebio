from django.db import migrations, models

def populate_compendia_results(apps, schema):
    """ Populate compendia_results with computed files """

    ComputedFile = apps.get_model('data_refinery_common', 'ComputedFile')
    CompendiaResult = apps.get_model('data_refinery_common', 'CompendiaResult')
    CompendiaResultOrganismAssociation = apps.get_model('data_refinery_common',
            'CompendiaResultOrganismAssociation')

    # get all existing compendia computed files
    compendia_computed_files = ComputedFile.objects.filter(is_compendia=True)

    # create as compendia result from each of the results
    compendia_results = []
    for computed_file in compendia_computed_files:
        compendia_result = CompendiaResult()
        compendia_result.quant_sf_only = computed_file.quant_sf_only
        compendia_result.compendia_version = computed_file.compenida_version
        compendia_result.svd_algorithm = computed_file.compendia_version
        compendia_result.result = computed_file.result
        compendia_result.primary_organism = computed_file.compendia_organism
        compendia_results.append(compendia_result)

    CompendiaResult.objects.bulk_create(compendia_results)

    compendia_result_organism_associations = []
    for compendia_result in compendia_results:
        compendia_result_organism_association = SampleComputedFileAssociation()
        compendia_result_organism_association.compendia_result = compendia_result
        compendia_result_organism_association.organism = compendia_result.primary_organism

    CompendiaResultOrganismAssociation.objects.bulk_create(
            compendia_result_organism_associations)

class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0042_auto_20191102_1641'),
    ]

    operations = [
        migrations.RunPython(populate_compendia_results)
    ]
