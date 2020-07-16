##
# Contains CreateAPITokenView, APITokenView, and needed serializer
##

from django.utils.decorators import method_decorator
from rest_framework import mixins, serializers, viewsets

from drf_yasg.utils import swagger_auto_schema

from data_refinery_common.models import APIToken


class APITokenSerializer(serializers.ModelSerializer):
    class Meta:
        model = APIToken
        fields = ("id", "is_activated", "terms_and_conditions")
        extra_kwargs = {
            "id": {"read_only": True},
            "is_activated": {"read_only": False},
            "terms_and_conditions": {"read_only": True},
        }


@method_decorator(
    name="create",
    decorator=swagger_auto_schema(
        operation_description="""
    This endpoint can be used to create and activate tokens. These tokens can be used
    in requests that provide urls to download computed files. They are a way to accept
    our terms of service.

    ```py
    import requests
    import json

    response = requests.post('https://api.refine.bio/v1/token/')
    token_id = response.json()['id']
    response = requests.put('https://api.refine.bio/v1/token/' + token_id + '/', json.dumps({'is_activated': True}), headers={'Content-Type': 'application/json'})
    ```

    The token id needs to be provided in the HTTP request in the API-KEY header.

    References
    - [https://github.com/AlexsLemonade/refinebio/issues/731]()
    - [https://github.com/AlexsLemonade/refinebio-frontend/issues/560]()
    """
    ),
)
@method_decorator(
    name="retrieve",
    decorator=swagger_auto_schema(operation_description="Return details about a specific token.",),
)
@method_decorator(
    name="update",
    decorator=swagger_auto_schema(
        operation_description="This can be used to activate a specific token by sending `is_activated: true`."
    ),
)
class APITokenView(
    mixins.CreateModelMixin,
    mixins.UpdateModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet,
):
    """
    Create, read, and modify Api Tokens.
    """

    model = APIToken
    lookup_field = "id"
    queryset = APIToken.objects.all()
    serializer_class = APITokenSerializer
