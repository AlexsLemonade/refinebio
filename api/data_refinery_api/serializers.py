from rest_framework import serializers
from data_refinery_common.models import (
    Experiment
)

class ExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ('id', 'title', 'description')