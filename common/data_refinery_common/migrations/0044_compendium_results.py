from django.db import migrations, models

def populate_compendium_results(apps, schema_editor):
    """ Fetch existing compendia and create compendium results"""
    CompendiumResult = apps.get_model('data_refinery_common', 'CompendiumResult')
    CompendiumResultOrganismAssociation = apps.get_model('data_refinery_common',
                                                         'CompendiumResultOrganismAssociation')
    ComputedFile = apps.get_modoel('data_refinery_common', 'ComputedFile')

    compendium_computed_files = ComputedFile.object.filter(is_compendia=True)

    compendium_results = []
    for computed_file in compendium_computed_files:
        compendium_result = CompendiumResult()
        compendium_result.quant_sf_only = computed_file.quant_sf_only
        compendium_result.svd_alorgithm = computed_file.svd_alorgithm
        compendium_result.primary_organism = computed_file.compendia_organism
        compendium_result.compendium_version = computed_file.compendia_version
        compendium_result.result = computed_file.result
        compendium_results.append(compendium_result)

    CompendiumResult.objects.bulk_create(compendium_results)

    # create compendium_result.organism relationships
    compendium_result_organism_associations = []

    for compendium_result in CompendiumResult.objects.all():
        compendium_result_organism_association = CompendiumResultOrganismAssociation()
        compendium_result_organism_association.compendium_result = compendium_result
        compendium_result_organism_association.organism = compendium_result.primary_organism
        compendium_result_organism_associations.append(
            compendium_result_organism_association)

    CompendiumResultOrganismAssociation.bulk_create(
            compendium_result_organism_associations)


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0043_auto_20191108_2144'),
    ]

    operations = [
        migrations.RunPython(set_computed_file_svd_algorithm_none),
    ]
