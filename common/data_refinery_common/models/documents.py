from elasticsearch_dsl import analyzer

from django_elasticsearch_dsl import DocType, Index, fields 
from elasticsearch_dsl.analysis import token_filter

from .models import Sample, Experiment, Organism

experiment_index = Index('experiments')
experiment_index.settings(
    number_of_shards=1,
    number_of_replicas=0,
    max_result_window=9999999
)

# via https://django-elasticsearch-dsl-drf.readthedocs.io/en/0.17.2/advanced_usage_examples.html?highlight=ngram#id8
# via https://github.com/barseghyanartur/django-elasticsearch-dsl-drf/issues/110
edge_ngram_completion_filter = token_filter(
    'edge_ngram_completion_filter',
    type="edge_ngram",
    min_gram=4,
    max_gram=12
)
html_strip = analyzer(
    'html_strip',
    tokenizer="standard",
    filter=["standard", "lowercase", "stop", "snowball", edge_ngram_completion_filter],
    char_filter=["html_strip"]
)
html_strip_no_ngram = analyzer(
    'html_strip_no_ngram',
    tokenizer="standard",
    filter=["standard", "lowercase", "stop", "snowball"],
    char_filter=["html_strip"]
)
standard = analyzer(
    'standard',
    tokenizer="keyword",
    filter=[],
)

@experiment_index.doc_type
class ExperimentDocument(DocType):
    """ Our Experiment ElasticSearch Document, which
    corresponds to our Experiment model. """

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
    technology = fields.TextField(
        analyzer=html_strip_no_ngram,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    organism_names = fields.TextField(
        analyzer=html_strip_no_ngram,
        fielddata=True,
        fields={'raw': fields.KeywordField()}
    )
    platform_names = fields.TextField(
        analyzer=standard,
        fielddata=True,
        fields={'raw': fields.TextField()}
    )
    platform_accession_codes = fields.TextField(
        analyzer=standard,
        fielddata=True,
        fields={'raw': fields.TextField()}
    )
    
    # Basic Fields
    accession_code = fields.TextField()
    alternate_accession_code = fields.TextField()
    submitter_institution = fields.TextField()
    publication_doi = fields.TextField()
    has_publication = fields.BooleanField()
    sample_metadata_fields = fields.TextField()
    pubmed_id = fields.TextField()
    num_total_samples = fields.IntegerField()
    num_processed_samples = fields.IntegerField()
    source_first_published = fields.DateField()

    # FK/M2M
    # We actually don't use any ForeignKeys in our Experiment document,
    # but if we did, we'd do it like this. The function `get_instances_from_related` is similarly required,
    # as is the `related_models` field in the Meta class.

    # organisms = fields.NestedField(properties={
    #     'name': fields.KeywordField(),
    #     'taxonomy_id': fields.IntegerField(),
    #     'pk': fields.IntegerField(),
    # })
    # 
    # def get_instances_from_related(self, related_instance):
    #     return related_instance.experts_set.all()

    class Meta:
        model = Experiment

        fields = [
           'id',
        ]
        # related_models = [Sample, Organism]
