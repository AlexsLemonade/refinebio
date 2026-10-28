##
# Experiment document view.
#
# Hand-rolled against opensearch-py + DRF: replaces the previous
# django-elasticsearch-dsl-drf DocumentViewSet + filter-backend stack that
# the OpenSearch migration abandoned.
##

from django.http import QueryDict
from django.utils.decorators import method_decorator
from rest_framework import serializers
from rest_framework.renderers import BrowsableAPIRenderer, JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView

from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import OpenApiParameter, extend_schema
from opensearchpy import A, Q

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_common.models.documents import ExperimentDocument

# ---------------------------------------------------------------------------
# Config — mirrors the previous DocumentViewSet class attributes exactly.
# Changing anything here changes the public API.
# ---------------------------------------------------------------------------

# {public search param: boost (or None for unboosted)}.
SEARCH_FIELDS = {
    "title": 10,
    "publication_authors": 8,  # "People will search themselves"
    "sample_keywords": 7,
    "publication_title": 5,
    "submitter_institution": 3,
    "description": 2,
    "accession_code": None,
    "alternate_accession_code": None,
    "publication_doi": None,
    "pubmed_id": None,
    "sample_metadata_fields": None,
    "platform_names": None,
}

# {public filter param: {"field": <es field>, "lookups": set(<suffixes>)}}.
# Suffixes supported: "" (exact term), "in" (comma-separated terms), "gt",
# "gte", "lt", "lte".
FILTER_FIELDS = {
    "id": {"field": "_id", "lookups": {"", "in"}},
    "technology": {"field": "technology", "lookups": {""}},
    "has_publication": {"field": "has_publication", "lookups": {""}},
    "accession_code": {"field": "accession_code", "lookups": {"", "in"}},
    "alternate_accession_code": {"field": "alternate_accession_code", "lookups": {"", "in"}},
    "platform": {"field": "platform_accession_codes", "lookups": {""}},
    "organism": {"field": "organism_names.raw", "lookups": {""}},
    "downloadable_organism": {"field": "downloadable_organism_names.raw", "lookups": {""}},
    "num_processed_samples": {
        "field": "num_processed_samples",
        "lookups": {"", "in", "gt", "gte", "lt", "lte"},
    },
    "num_downloadable_samples": {
        "field": "num_downloadable_samples",
        "lookups": {"", "in", "gt", "gte", "lt", "lte"},
    },
    "sample_keywords": {"field": "sample_keywords", "lookups": {""}},
}

# {public ordering alias: es field}.
ORDERING_FIELDS = {
    "id": "id",
    "title": "title.raw",
    "description": "description.raw",
    "num_total_samples": "num_total_samples",
    "num_downloadable_samples": "num_downloadable_samples",
    "source_first_published": "source_first_published",
}

DEFAULT_ORDERING = (
    "_score",
    "-num_total_samples",
    "id",
    "title",
    "description",
    "-source_first_published",
)

# Facets kept as (public key, es field, options) tuples. All four are
# enabled by default; "has_publication" uses the key_as_string shape in the
# transform below.
FACETS = (
    ("technology", "technology", {}),
    ("downloadable_organism_names", "downloadable_organism_names.raw", {"size": 999999}),
    ("platform_accession_codes", "platform_accession_codes", {"size": 999999}),
    ("has_publication", "has_publication", {}),
)

CARDINALITY_FIELD = "downloadable_samples"
CARDINALITY_METRIC = "total_samples"
CARDINALITY_PRECISION = 40000

DEFAULT_LIMIT = 10
MAX_LIMIT = 10000


# ---------------------------------------------------------------------------
# Serializer — flat DRF Serializer, replaces the old DocumentSerializer.
# ---------------------------------------------------------------------------


class ExperimentDocumentSerializer(serializers.Serializer):
    """Serializes a hit from the Experiment index."""

    id = serializers.CharField(read_only=True)
    title = serializers.CharField(read_only=True)
    publication_title = serializers.CharField(read_only=True)
    description = serializers.CharField(read_only=True)
    technology = serializers.CharField(read_only=True)
    accession_code = serializers.CharField(read_only=True)
    alternate_accession_code = serializers.CharField(read_only=True)
    submitter_institution = serializers.CharField(read_only=True)
    has_publication = serializers.BooleanField(read_only=True)
    publication_doi = serializers.CharField(read_only=True)
    publication_authors = serializers.ListField(child=serializers.CharField(), read_only=True)
    sample_metadata_fields = serializers.ListField(child=serializers.CharField(), read_only=True)
    platform_names = serializers.ListField(child=serializers.CharField(), read_only=True)
    platform_accession_codes = serializers.ListField(child=serializers.CharField(), read_only=True)
    organism_names = serializers.ListField(child=serializers.CharField(), read_only=True)
    downloadable_organism_names = serializers.ListField(
        child=serializers.CharField(), read_only=True
    )
    pubmed_id = serializers.CharField(read_only=True)
    num_total_samples = serializers.IntegerField(read_only=True)
    num_processed_samples = serializers.IntegerField(read_only=True)
    num_downloadable_samples = serializers.IntegerField(read_only=True)
    source_first_published = serializers.DateTimeField(read_only=True)


# ---------------------------------------------------------------------------
# Renderer — the old view suppressed form rendering so the browsable API
# would not try to render a POST form it cannot build. Keep it.
# ---------------------------------------------------------------------------


class FormlessBrowsableAPIRenderer(BrowsableAPIRenderer):
    def show_form_for_method(self, view, method, request, instance):
        return False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _coerce_param_value(value):
    """Normalise a POST body value into a list of strings, matching the old
    POSTFilteringFilterBackend.MockRequest adapter."""
    if isinstance(value, list):
        return [_coerce_param_value(v)[0] for v in value]
    if isinstance(value, str):
        return [value]
    if isinstance(value, (int, bool)):
        return [str(value)]
    raise InvalidFilters(
        message="Invalid type {} for filter value {}".format(type(value), str(value))
    )


def _request_params(request):
    """Return a QueryDict holding the filter/search params, pulled from the
    query string for GET and the JSON body for POST."""
    if request.method == "POST":
        params = QueryDict(mutable=True)
        for key, value in request.data.items():
            for v in _coerce_param_value(value):
                params.appendlist(key, v)
        return params
    return request.query_params


def _split_filter_key(key):
    """Return (public_name, lookup_suffix) for a filter query-param key."""
    if "__" in key:
        name, _, suffix = key.rpartition("__")
        return name, suffix
    return key, ""


def _build_filter_query(field, suffix, values):
    """Build the opensearch-py Q for one filter param."""
    if suffix == "in":
        terms = []
        for raw in values:
            terms.extend(v for v in raw.split(",") if v)
        return Q("terms", **{field: terms})
    if suffix in {"gt", "gte", "lt", "lte"}:
        return Q("range", **{field: {suffix: values[0]}})
    # Exact term(s).
    if len(values) == 1:
        return Q("term", **{field: values[0]})
    return Q("terms", **{field: list(values)})


def _apply_filters(search, params):
    """Translate the public filter params into opensearch-py .filter() calls."""
    for key in params:
        name, suffix = _split_filter_key(key)
        config = FILTER_FIELDS.get(name)
        if not config:
            continue
        if suffix not in config["lookups"]:
            continue
        values = params.getlist(key)
        search = search.filter(_build_filter_query(config["field"], suffix, values))
    return search


def _apply_free_text_search(search, params):
    """Apply the free-text `search=` param as a multi_match across the boosted
    field set."""
    text = params.get("search")
    if not text:
        return search
    fields = []
    for field, boost in SEARCH_FIELDS.items():
        fields.append(f"{field}^{boost}" if boost else field)
    return search.query("multi_match", query=text, fields=fields)


def _apply_ordering(search, params):
    """Apply `ordering=field` (one or many, prefix `-` for desc), falling back
    to the default when none is given."""
    orderings = params.getlist("ordering") if hasattr(params, "getlist") else []
    if orderings:
        resolved = []
        for raw in orderings:
            desc = raw.startswith("-")
            name = raw.lstrip("-")
            field = ORDERING_FIELDS.get(name)
            if field:
                resolved.append(f"-{field}" if desc else field)
        if resolved:
            return search.sort(*resolved)
    return search.sort(*DEFAULT_ORDERING)


def _apply_facets(search):
    """Attach the four default facets, each with the custom cardinality
    sub-metric on downloadable_samples."""
    for key, field, options in FACETS:
        agg = A("terms", field=field, **options)
        agg.metric(
            CARDINALITY_METRIC,
            "cardinality",
            field=CARDINALITY_FIELD,
            precision_threshold=CARDINALITY_PRECISION,
        )
        search.aggs.bucket(key, agg)
    return search


def _transform_facets(aggregations):
    """Reshape raw aggregation buckets into the API's facet envelope:
    {field: {bucket_key: total_samples}}."""
    result = {}
    for key, _field, _options in FACETS:
        agg = aggregations.get(key, {})
        buckets = agg.get("buckets", [])
        group = {}
        for bucket in buckets:
            # has_publication is a boolean term; use the stringified key so
            # consumers can address "true"/"false" directly.
            label = bucket.get("key_as_string") if key == "has_publication" else bucket.get("key")
            group[label] = bucket[CARDINALITY_METRIC]["value"]
        result[key] = group
    return result


def _hit_to_dict(hit):
    """Convert a Hit into the dict the serializer expects. `meta.id` surfaces
    as the top-level `id` the old view relied on."""
    data = hit.to_dict()
    data.setdefault("id", hit.meta.id)
    return data


def _parse_pagination(params):
    """Return (limit, offset) with the same defaults the old
    ESLimitOffsetPagination used."""
    try:
        limit = int(params.get("limit", DEFAULT_LIMIT))
    except (TypeError, ValueError):
        limit = DEFAULT_LIMIT
    try:
        offset = int(params.get("offset", 0))
    except (TypeError, ValueError):
        offset = 0
    limit = max(0, min(limit, MAX_LIMIT))
    offset = max(0, offset)
    return limit, offset


# ---------------------------------------------------------------------------
# View
# ---------------------------------------------------------------------------


@method_decorator(
    name="get",
    decorator=extend_schema(
        parameters=[
            OpenApiParameter(
                name="accession_code",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.STR,
                description="Filter by accession code; may be repeated.",
            ),
            OpenApiParameter(
                name="technology",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.STR,
                description="Filter by technology; may be repeated.",
            ),
            OpenApiParameter(
                name="has_publication",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.STR,
                description="Filter the results that have associated publications with `?has_publication=true`.",
            ),
            OpenApiParameter(
                name="platform",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.STR,
                description="Filter by platform; may be repeated.",
            ),
            OpenApiParameter(
                name="organism",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.STR,
                description="Filter by organism; may be repeated.",
            ),
            OpenApiParameter(
                name="num_processed_samples",
                location=OpenApiParameter.QUERY,
                type=OpenApiTypes.NUMBER,
                description="Filter by number of processed samples; supports __gt, __gte, __lt, __lte.",
            ),
        ],
        description="""
Use this endpoint to search among the experiments.

Powered by OpenSearch. Filter parameters support the suffixes `__in`, `__gt`,
`__gte`, `__lt`, and `__lte`.

The response includes a `facets` field with per-filter bucket counts.

Example Requests:
```
?search=medulloblastoma
?id=1
?search=medulloblastoma&technology=microarray&has_publication=true
?ordering=source_first_published
```

This endpoint also accepts POST requests for larger queries. Any of the
filters accepted as query parameters are also accepted in a JSON object in
the request body.

Example Requests (from our tests):
```python
import requests
import json

headers = {
    'Content-Type': 'application/json',
}

# Basic filter
search = {"accession_code": "GSE123"}
requests.post(host + '/v1/search/', json.dumps(search), headers=headers)

# __in filter
search = {"accession_code__in": ["GSE123"]}
requests.post(host + '/v1/search/', json.dumps(search), headers=headers)

# numeric filter
search = {"num_downloadable_samples__gt": 0}
requests.post(host + '/v1/search/', json.dumps(search), headers=headers)
```
""",
    ),
)
class ExperimentDocumentView(APIView):
    """OpenSearch-backed experiment search."""

    renderer_classes = [JSONRenderer, FormlessBrowsableAPIRenderer]
    serializer_class = ExperimentDocumentSerializer

    def get(self, request, *args, **kwargs):
        return self._search(request)

    @extend_schema(exclude=True)
    def post(self, request, *args, **kwargs):
        return self._search(request)

    def _search(self, request):
        params = _request_params(request)
        limit, offset = _parse_pagination(params)

        search = ExperimentDocument.search()
        search = _apply_free_text_search(search, params)
        search = _apply_filters(search, params)
        search = _apply_ordering(search, params)
        _apply_facets(search)

        # Force eager total so hits.total.value is populated on OS 7+/3.x
        # defaults.
        search = search.extra(track_total_hits=True)

        response = search[offset : offset + limit].execute()
        aggregations = response.aggregations.to_dict() if response.aggregations else {}

        serializer = self.serializer_class([_hit_to_dict(hit) for hit in response.hits], many=True)

        total = response.hits.total.value
        return Response(
            {
                "count": total,
                "next": self._page_url(request, limit, offset + limit)
                if offset + limit < total
                else None,
                "previous": self._page_url(request, limit, max(0, offset - limit))
                if offset > 0
                else None,
                "results": serializer.data,
                "facets": _transform_facets(aggregations),
            }
        )

    @staticmethod
    def _page_url(request, limit, offset):
        """Build a limit/offset-swapped URL for the next/previous links."""
        query = request.query_params.copy()
        query["limit"] = str(limit)
        query["offset"] = str(offset)
        return request.build_absolute_uri(f"{request.path}?{query.urlencode()}")
