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
    ComputedFile
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
                    'accession_code', 
                    'organism', 
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
                    'accession_code',
                    'organism',
                    'annotations',
                    'results',
                    'source_archive_url',
                    'has_raw',
                    'is_downloaded',
                    'is_processed',
                    'created_at',
                    'last_modified',
                )

##
# Experiments
##

class ExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = (
                    'id',
                    'title',
                    'description',
                    'accession_code',
                    'platform_accession_code',
                    'source_database',
                    'source_url',
                    'platform_name',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'pubmed_id',
                    'samples',
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
                    'platform_name',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'pubmed_id',
                    'source_first_published',
                    'source_last_updated',
                    'platform_accession_code',
                    'submitter_institution',
                    'last_modified',
                    'created_at',
                )

class PlatformSerializer(serializers.ModelSerializer):

    class Meta:
        model = Experiment
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
                    'has_raw',
                    'is_downloaded',
                    'is_processed'
                )

##
# Jobs
##

class SurveyJobSerializer(serializers.ModelSerializer):

    class Meta:
        model = SurveyJob
        fields = (
                    'id',
                    'num_retries',
                    'retried',
                    'worker_id',
                    'worker_version',
                    'failure_reason',
                    'source_type',
                    'success',
                    'start_time',
                    'end_time'
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
                    'end_time'
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
                    'end_time'
                )

