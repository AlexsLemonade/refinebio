##
# Contains the views QNTargetsAvailable and QNTargetsDetailView
##

from django.utils.decorators import method_decorator
from rest_framework import generics, serializers
from rest_framework.exceptions import NotFound

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_api.views.relation_serializers import (
    ComputationalResultNoFilesRelationSerializer,
    OrganismRelationSerializer,
)
from data_refinery_common.models import ComputationalResultAnnotation, ComputedFile, Organism


class QNTargetSerializer(serializers.ModelSerializer):
    result = ComputationalResultNoFilesRelationSerializer(many=False)

    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "is_qn_target",
            "sha1",
            "s3_bucket",
            "s3_key",
            "s3_url",
            "created_at",
            "last_modified",
            "result",
        )


class QNTargetsAvailable(generics.ListAPIView):
    """
    This is a list of all of the organisms which have available QN Targets
    """

    serializer_class = OrganismRelationSerializer
    paginator = None

    def get_queryset(self):
        return Organism.get_objects_with_qn_targets()


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="organism_name",
                in_=openapi.IN_PATH,
                type=openapi.TYPE_STRING,
                description="Eg `DANIO_RERIO`, `MUS_MUSCULUS`",
            )
        ],
        responses={404: "QN Target not found for the given organism."},
    ),
)
class QNTargetsDetailView(generics.RetrieveAPIView):
    """
    Get a detailed view of the Quantile Normalization file for an organism.
    """

    serializer_class = QNTargetSerializer

    def get_object(self):
        organism = self.kwargs["organism_name"]
        organism = organism.upper().replace(" ", "_")
        try:
            organism_id = Organism.get_object_for_name(organism).id
            annotation = (
                ComputationalResultAnnotation.objects.filter(
                    data__organism_id=organism_id, data__is_qn=True
                )
                .order_by("-created_at")
                .first()
            )
            qn_target = annotation.result.computedfile_set.first()
        except Exception:
            raise NotFound("Don't have a target for that organism!")
        if not qn_target:
            raise NotFound("Don't have a target for that organism!!")
        return qn_target
