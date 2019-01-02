from elasticsearch_dsl import analyzer

from django_elasticsearch_dsl import DocType, Index, fields 

from .models import Sample, Experiment

experiment_index = Index('experiments')
experiment_index.settings(
    number_of_shards=1,
    number_of_replicas=0
)

html_strip = analyzer(
    'html_strip',
    tokenizer="standard",
    filter=["standard", "lowercase", "stop", "snowball"],
    char_filter=["html_strip"]
)

@experiment_index.doc_type
class ExperimentDocument(DocType):
    """Experiment elasticsearch document"""

    samples = fields.NestedField(properties={
        'title': fields.StringField()
    })

    class Meta:
        model = Experiment

        fields = [
            'id',
            'title',
            'description'
        ]
