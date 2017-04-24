from enum import Enum
from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel
from data_refinery_models.models.surveys import SurveyJob


class BatchStatuses(Enum):
    NEW = "NEW"
    DOWNLOADED = "DOWNLOADED"
    PROCESSED = "PROCESSED"


class Batch(TimeTrackedModel):
    survey_job = models.ForeignKey(SurveyJob)
    source_type = models.CharField(max_length=256)
    size_in_bytes = models.IntegerField()
    download_url = models.CharField(max_length=2048)
    raw_format = models.CharField(max_length=256, null=True)
    processed_format = models.CharField(max_length=256, null=True)
    pipeline_required = models.CharField(max_length=256)
    platform_accession_code = models.CharField(max_length=32)
    experiment_accession_code = models.CharField(max_length=32)
    experiment_title = models.CharField(max_length=32)
    status = models.CharField(max_length=20)
    release_date = models.DateTimeField()
    last_uploaded_date = models.DateTimeField()
    # api Revision?!? -- if so probably better as a KV record
    file_name = models.CharField(max_length=1024)

    # This field will denote where in our system the file can be found
    internal_location = models.CharField(max_length=256, null=True)

    # This corresponds to the organism taxonomy ID from NCBI
    # This looks like it may not always be easy... but
    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=homo+sapiens&field=scin # noqa
    # will get us the ID.
    # This is probably worth more discussion. However an idea we could play
    # with is creating some kind of internal database that is caching
    # organism code -> scientific name
    # This would allow us to beef up the discovery of organism codes if there
    # are organisms which aren't returned by that query.
    # If our code cannot figure it out then it could store the organism_name,
    # set the organism_id to 0, and then notify us so we can figure out how
    # to handle the new special case.
    # Initial data model idea looks like:
    #     organism_id -- NCBI taxonomy ID
    #     organism_name -- name of the organism used to lookup the id
    #     scientific_name -- boolean indicating if this is the scientific name
    organism_id = models.IntegerField()
    organism_name = models.CharField(max_length=256)

    class Meta:
        db_table = "batches"


class BatchKeyValue(TimeTrackedModel):
    """
    This table is used for tracking fields onto a Batch record that would
    be sparsely populated if it was its own column.
    I.e. one source may have an extra field or two that are worth tracking
    but are specific to that source.
    """
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "batch_key_values"
