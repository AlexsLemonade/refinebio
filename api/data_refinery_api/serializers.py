from rest_framework import serializers
from rest_framework_hstore.fields import HStoreField
from data_refinery_common.models import ProcessorJob, DownloaderJob, SurveyJob
from data_refinery_common.models import (
    Experiment,
    ExperimentAnnotation,
    Sample,
    SampleAnnotation,
    Organism,
    OriginalFile,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset
)

##
# Organism
##

class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = (
                    'name',
                    'taxonomy_id',
                )

##
# Results
##

class ComputationalResultAnnotationSerializer(serializers.ModelSerializer):
    data = HStoreField()

    class Meta:
        model = ComputationalResultAnnotation
        fields = (
                    'id',
                    'data',
                    'is_ccdl',
                    'created_at',
                    'last_modified'
                )

class ComputedFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputedFile
        fields = (
                    'id',
                    'filename',
                    'size_in_bytes',
                    'sha1',
                    's3_bucket',
                    's3_key',
                    'created_at',
                    'last_modified'
                )

class ComputationalResultSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationSerializer(many=True, source='computationalresultannotation_set')
    files = ComputedFileSerializer(many=True, source='computedfile_set')

    class Meta:
        model = ComputationalResult
        fields = (
                    'id',
                    'command_executed',
                    'program_version',
                    'system_version',
                    'is_ccdl',
                    'annotations',
                    'files',
                    'time_start',
                    'time_end',
                    'created_at',
                    'last_modified'
                )

##
# Samples
##

class SampleSerializer(serializers.ModelSerializer):
    organism = OrganismSerializer(many=False)

    class Meta:
        model = Sample
        fields = (
                    'id',
                    'title',
                    'accession_code',
                    'source_database',
                    'organism',
                    'platform_accession_code',
                    'platform_name',
                    'technology',
                    'manufacturer',
                    'is_downloaded',
                    'is_processed',
                    'created_at',
                    'last_modified',
                )

class SampleAnnotationSerializer(serializers.ModelSerializer):
    data = HStoreField()

    class Meta:
        model = SampleAnnotation
        fields = (
                    'data',
                    'is_ccdl',
                    'created_at',
                    'last_modified',
                )

class DetailedSampleSerializer(serializers.ModelSerializer):
    annotations = SampleAnnotationSerializer(many=True, source='sampleannotation_set')
    organism = OrganismSerializer(many=False)
    results = ComputationalResultSerializer(many=True)

    class Meta:
        model = Sample
        fields = (
                    'id',
                    'title',
                    'accession_code',
                    'source_database',
                    'organism',
                    'platform_accession_code',
                    'platform_name',
                    'technology',
                    'manufacturer',
                    'annotations',
                    'results',
                    'pipelines',
                    'source_archive_url',
                    'has_raw',
                    'sex',
                    'age',
                    'specimen_part',
                    'genotype',
                    'disease',
                    'disease_stage',
                    'cell_line',
                    'treatment',
                    'race',
                    'subject',
                    'compound',
                    'time',
                    'is_downloaded',
                    'is_processed',
                    'created_at',
                    'last_modified',
                )

##
# Experiments
##

class ExperimentSerializer(serializers.ModelSerializer):
    organisms = serializers.StringRelatedField(many=True)
    platforms = serializers.ReadOnlyField()
    samples = serializers.StringRelatedField(many=True)

    class Meta:
        model = Experiment
        fields = (
                    'id',
                    'title',
                    'description',
                    'accession_code',
                    'source_database',
                    'source_url',
                    'platforms',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'pubmed_id',
                    'samples',
                    'organisms',
                    'submitter_institution',
                    'created_at',
                    'last_modified'
                )

class ExperimentAnnotationSerializer(serializers.ModelSerializer):
    data = HStoreField()

    class Meta:
        model = ExperimentAnnotation
        fields = (
                    'data',
                    'is_ccdl',
                    'created_at',
                    'last_modified',
                )

class DetailedExperimentSerializer(serializers.ModelSerializer):
    annotations = ExperimentAnnotationSerializer(many=True, source='experimentannotation_set')
    samples = SampleSerializer(many=True)

    class Meta:
        model = Experiment
        fields = (
                    'id',
                    'title',
                    'description',
                    'annotations',
                    'samples',
                    'protocol_description',
                    'accession_code',
                    'source_database',
                    'source_url',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'pubmed_id',
                    'source_first_published',
                    'source_last_modified',
                    'submitter_institution',
                    'last_modified',
                    'created_at',
                )

class PlatformSerializer(serializers.ModelSerializer):

    class Meta:
        model = Sample
        fields = (
                    'platform_accession_code',
                    'platform_name',
                )

class InstitutionSerializer(serializers.ModelSerializer):

    class Meta:
        model = Experiment
        fields = (
                    'submitter_institution',
                )

##
# Files
##

class OriginalFileSerializer(serializers.ModelSerializer):

    class Meta:
        model = OriginalFile
        fields = (
                    'id',
                    'filename',
                    'size_in_bytes',
                    'sha1',
                    'source_url',
                    'source_filename',
                    'is_downloaded',
                    'is_archive',
                    'has_raw',
                    'is_downloaded',
                    'created_at',
                    'last_modified'
                )

##
# Jobs
##

class SurveyJobSerializer(serializers.ModelSerializer):

    class Meta:
        model = SurveyJob
        fields = (
                    'id',
                    'source_type',
                    'success',
                    'start_time',
                    'end_time',
                    'created_at',
                    'last_modified'
                )

class DownloaderJobSerializer(serializers.ModelSerializer):
    original_files = OriginalFileSerializer(many=True)

    class Meta:
        model = DownloaderJob
        fields = (
                    'id',
                    'downloader_task',
                    'num_retries',
                    'retried',
                    'worker_id',
                    'worker_version',
                    'failure_reason',
                    'success',
                    'original_files',
                    'start_time',
                    'end_time',
                    'created_at',
                    'last_modified'
                )

class ProcessorJobSerializer(serializers.ModelSerializer):
    original_files = OriginalFileSerializer(many=True)

    class Meta:
        model = ProcessorJob
        fields = (
                    'id',
                    'pipeline_applied',
                    'num_retries',
                    'retried',
                    'worker_id',
                    'worker_version',
                    'failure_reason',
                    'success',
                    'original_files',
                    'start_time',
                    'end_time',
                    'created_at',
                    'last_modified'
                )

##
# Datasets
##

def validate_dataset(data):
    """ Basic dataset validation. Currently only checks formatting, not values. """
    if data['data'] != None:
        if type(data['data']) != dict:
            raise serializers.ValidationError("`data` must be a dict of lists.")

        for key, value in data['data'].items():
            if type(value) != list:
                raise serializers.ValidationError("`data` must be a dict of lists. Problem with `" + str(key) + "`")

    else:
        raise serializers.ValidationError("`data` must be a dict of lists.")

class CreateDatasetSerializer(serializers.ModelSerializer):

    class Meta:
        model = Dataset
        fields = (
                    'id',
                    'data',
                    'email_address'
            )

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        try:
            validate_dataset(data)
        except Exception:
            raise
        return data

class DatasetSerializer(serializers.ModelSerializer):

    start = serializers.NullBooleanField(required=False)

    class Meta:
        model = Dataset
        fields = (
                    'id',
                    'data',
                    'aggregate_by',
                    'is_processing',
                    'is_processed',
                    'is_available',
                    'email_address',
                    'expires_on',
                    's3_bucket',
                    's3_key',
                    'created_at',
                    'last_modified',
                    'start'
            )
        extra_kwargs = {
                        'id': {
                            'read_only': True,
                        },
                        'is_processing': {
                            'read_only': True,
                        },
                        'is_processed': {
                            'read_only': True,
                        },
                        'is_available': {
                            'read_only': True,
                        },
                        'expires_on': {
                            'read_only': True,
                        },
                        's3_bucket': {
                            'read_only': True,
                        },
                        's3_key': {
                            'read_only': True,
                        },
                        'created_at': {
                            'read_only': True,
                        },
                        'last_modified': {
                            'read_only': True,
                        }
                    }

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        try:
            validate_dataset(data)
        except Exception:
            raise
        return data
