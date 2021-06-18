from django.test import TestCase, tag


class ImportTestCase(TestCase):
    def test_imports(self):
        # Make sure we can import the foreman tests
        import data_refinery_foreman.foreman.management.commands.test_assoc_experiment_results
        import data_refinery_foreman.foreman.management.commands.test_correct_affy_cdfs
        import data_refinery_foreman.foreman.management.commands.test_create_compendia
        import data_refinery_foreman.foreman.management.commands.test_create_missing_downloader_jobs
        import data_refinery_foreman.foreman.management.commands.test_create_missing_processor_jobs
        import data_refinery_foreman.foreman.management.commands.test_create_quantpendia
        import data_refinery_foreman.foreman.management.commands.test_import_external_sample_attributes
        import data_refinery_foreman.foreman.management.commands.test_import_external_sample_keywords
        import data_refinery_foreman.foreman.management.commands.test_organism_shepherd
        import data_refinery_foreman.foreman.management.commands.test_rerun_salmon_old_samples
        import data_refinery_foreman.foreman.management.commands.test_retry_samples
        import data_refinery_foreman.foreman.management.commands.test_run_tximport
        import data_refinery_foreman.foreman.management.commands.test_update_experiment_metadata
        import data_refinery_foreman.foreman.test_downloader_job_manager
        import data_refinery_foreman.foreman.test_end_to_end
        import data_refinery_foreman.foreman.test_job_control
        import data_refinery_foreman.foreman.test_job_requeuing
        import data_refinery_foreman.foreman.test_processor_job_manager
        import data_refinery_foreman.foreman.test_survey_job_manager
        import data_refinery_foreman.foreman.test_utils
        import data_refinery_foreman.surveyor.management.commands.test_unsurvey
        import data_refinery_foreman.surveyor.test_array_express
        import data_refinery_foreman.surveyor.test_external_source
        import data_refinery_foreman.surveyor.test_geo
        import data_refinery_foreman.surveyor.test_harmony
        import data_refinery_foreman.surveyor.test_sra
        import data_refinery_foreman.surveyor.test_surveyor
        import data_refinery_foreman.surveyor.test_transcriptome_index
