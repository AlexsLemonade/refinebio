from rest_framework import serializers
from data_refinery_common.models import (
    Experiment
)

class ExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = (	
        			'id', 
        			'title', 
        			'description',
        			'platform_accession_code',
        			'submitter_institution',
      				'created_at',
      				'last_modified',
        		)

class DetailedExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = (	
        			'id', 
        			'title', 
        			'description',
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