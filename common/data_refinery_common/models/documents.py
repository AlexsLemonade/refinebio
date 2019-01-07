from elasticsearch_dsl import analyzer

from django_elasticsearch_dsl import DocType, Index, fields 

from .models import Sample, Experiment, Organism

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

    # Keyword Fields
    title = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    publication_title = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    description = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    publication_authors = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )

    # Basic Fields
    technology = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    organism_names = fields.TextField(
        analyzer=html_strip,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )

    accession_code = fields.TextField()
    alternate_accession_code = fields.TextField()
    submitter_institution = fields.TextField()
    publication_doi = fields.TextField()
    has_publication = fields.BooleanField()
    sample_metadata_fields = fields.TextField()
    pubmed_id = fields.TextField()
    num_total_samples = fields.IntegerField()
    num_processed_samples = fields.IntegerField()

    # FK/M2M
    # organisms = fields.NestedField(properties={
    #     'name': fields.KeywordField(),
    #     'taxonomy_id': fields.IntegerField(),
    #     'pk': fields.IntegerField(),
    # })

    class Meta:
        model = Experiment

        fields = [
           'id',
        ]
        related_models = [Sample, Organism]
