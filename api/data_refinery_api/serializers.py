from rest_framework import serializers
from rest_framework_hstore.fields import HStoreField
from data_refinery_common.models import ProcessorJob, DownloaderJob, SurveyJob
from data_refinery_common.models import (
    Experiment,
    Sample,
    SampleAnnotation, 
    Organism,
    OriginalFile,
    ComputationalResult
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

class ComputationalResultSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputationalResult
        fields = (    
                    'id', 
                    'command_executed',
                    'program_version',
                    'system_version',
                    'is_ccdl',
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
                    'is_downloaded',
                    'is_processed',
                    'created_at',
                    'last_modified',
                    'annotations',
                    'results'
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
                    'platform_accession_code',
                    'samples',
                    'submitter_institution',
                    'created_at',
                    'last_modified'
                )

class DetailedExperimentSerializer(serializers.ModelSerializer):
    samples = SampleSerializer(many=True)

    class Meta:
        model = Experiment
        fields = (    
                    'id', 
                    'title', 
                    'description',
                    'samples',
                    'protocol_description',
                    'platform_accession_code',
                    'platform_name',
                    'has_publication',
                    'publication_title',
                    'publication_doi',
                    'pubmed_id',
                    'source_first_published',
                    'source_last_updated',
                    'platform_accession_code',
                    'submitter_institution',
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
                    'absolute_file_path',
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
                    'success',
                    'failure_reason',
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
                    'success',
                    'failure_reason',
                    'original_files',
                    'start_time',
                    'end_time'
                )

