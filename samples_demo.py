from django.db import connection
from data_refinery_common.models import (Experiment, Sample, ComputationalResult)
from django.db.models import Prefetch


experiment_id = Experiment.objects.values("id").get(accession_code="GSE11611")["id"]

base_queryset = (
    Sample.public_objects.prefetch_related("organism")
    .prefetch_related(
        Prefetch("results", queryset=ComputationalResult.objects.order_by("time_start"))
    )
    .prefetch_related("results__processor")
    .prefetch_related("results__computationalresultannotation_set")
    .prefetch_related("results__computedfile_set")
    .filter(experiments__in=[experiment_id])
)
