from rest_framework import serializers
from data_refinery_common.models import ProcessorJob, DownloaderJob, SurveyJob
from data_refinery_common.models import (
    Experiment,
    ExperimentAnnotation,
    Sample,
    SampleAnnotation,
    Organism,
    OrganismIndex,
    OriginalFile,
    Processor,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    APIToken
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
# Processor
##

class ProcessorSerializer(serializers.ModelSerializer):
    class Meta:
        model = Processor
        fields = (
            'id',
            'name',
            'version',
            'docker_image',
            'environment'
        )


##
# Transcriptome Index
##

class OrganismIndexSerializer(serializers.ModelSerializer):
    class Meta:
        model = OrganismIndex
        fields = (
                    's3_url',
                    'source_version',
                    'assembly_name',
                    'salmon_version',
                    'last_modified',
                )

##
# Results
##

class ComputationalResultAnnotationSerializer(serializers.ModelSerializer):

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
                    'is_smashable',
                    'is_qc',
                    'sha1',
                    's3_bucket',
                    's3_key',
                    'created_at',
                    'last_modified'
                )

class ComputationalResultSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationSerializer(many=True, source='computationalresultannotation_set')
    files = ComputedFileSerializer(many=True, source='computedfile_set')
    processor = ProcessorSerializer(many=False)

    class Meta:
        model = ComputationalResult
        fields = (
                    'id',
                    'commands',
                    'processor',
                    'is_ccdl',
                    'annotations',
                    'files',
                    'organism_index',
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
                    'pretty_platform',
                    'technology',
                    'manufacturer',
                    'protocol_info',
                    'is_processed',
                    'created_at',
                    'last_modified',
                )

class SampleAnnotationSerializer(serializers.ModelSerializer):

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
                    'pretty_platform',
                    'technology',
                    'manufacturer',
                    'protocol_info',
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
    processed_samples = serializers.StringRelatedField(many=True)
    pretty_platforms = serializers.ReadOnlyField()
    sample_metadata = serializers.ReadOnlyField(source='get_sample_metadata_fields')

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
                    'pretty_platforms',
                    'processed_samples',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'publication_authors',
                    'pubmed_id',
                    'samples',
                    'organisms',
                    'submitter_institution',
                    'created_at',
                    'last_modified',
                    'sample_metadata',
                )

class ExperimentAnnotationSerializer(serializers.ModelSerializer):

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
    organisms = OrganismSerializer(many=True)
    sample_metadata = serializers.ReadOnlyField(source='get_sample_metadata_fields')

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
                    'publication_authors',
                    'pubmed_id',
                    'source_first_published',
                    'source_last_modified',
                    'submitter_institution',
                    'last_modified',
                    'created_at',
                    'organisms',
                    'sample_metadata',
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
                    'ram_amount',
                    'volume_index',
                    'worker_version',
                    'failure_reason',
                    'success',
                    'original_files',
                    'datasets',
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

            try:
                if len(value) != len(set(value)):
                    raise serializers.ValidationError("Duplicate values detected in " + str(value))
            except Exception as e:
                raise serializers.ValidationError("Received bad dataset data: " + str(e))

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
                    'scale_by',
                    'is_processing',
                    'is_processed',
                    'is_available',
                    'has_email',
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

class DatasetDetailsSampleSerializer(serializers.ModelSerializer):
    """ This serializer contains all of the information about a sample needed for the download page
    """
    organism = OrganismSerializer(many=False)

    class Meta:
        model = Sample
        fields = (
                    'accession_code',
                    'organism',
                )

class DatasetDetailsExperimentSerializer(serializers.ModelSerializer):
    """ This serializer contains all of the information about an experiment needed for the download
    page
    """
    organisms = serializers.StringRelatedField(many=True)
    samples = serializers.StringRelatedField(many=True)
    sample_metadata = serializers.ReadOnlyField(source='get_sample_metadata_fields')
    pretty_platforms = serializers.ReadOnlyField()

    class Meta:
        model = Experiment
        fields = (
                    'id',
                    'title',
                    'accession_code',
                    'pretty_platforms',
                    'samples',
                    'organisms',
                    'sample_metadata',
                )

class DatasetDetailsSerializer(serializers.ModelSerializer):
    """ This serializer contains all of the information about a dataset needed for the download page
    """
    samples = DatasetDetailsSampleSerializer(read_only=True, many=True, source='get_samples')
    experiments = DatasetDetailsExperimentSerializer(read_only=True, many=True, source='get_experiments')

    class Meta:
        model = Dataset
        fields = (
                    'data',
                    'aggregate_by',
                    'scale_by',
                    'is_processing',
                    'is_processed',
                    'experiments',
                    'samples'
            )
        extra_kwargs = {
                        'is_processing': {
                            'read_only': True,
                        },
                        'is_processed': {
                            'read_only': True,
                        },
                    }

class APITokenSerializer(serializers.ModelSerializer):

    class Meta:
        model = APIToken
        fields = (
                'id',
                'is_activated',
                'terms_and_conditions'
            )
        extra_kwargs = {
                        'id': {
                            'read_only': True
                        },
                        'is_activated': {
                            'read_only': False
                        },
                        'terms_and_conditions': {
                            'read_only': True
                        }
                    }
