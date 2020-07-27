##
# Contains CreateAPITokenView, APITokenView, and needed serializer
##

from django.utils.decorators import method_decorator
from rest_framework import generics, serializers

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


class CreateAPITokenView(generics.CreateAPIView):
    """
    token_create

    This endpoint can be used to create and activate tokens. These tokens can be used
    in requests that provide urls to download computed files. Setting `is_activated` to
    true indicates agreement with refine.bio's [Terms of Use](https://www.refine.bio/terms)
    and [Privacy Policy](https://www.refine.bio/privacy).

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

    model = APIToken
    serializer_class = APITokenSerializer


@method_decorator(name="patch", decorator=swagger_auto_schema(auto_schema=None))
class APITokenView(generics.RetrieveUpdateAPIView):
    """
    Read and modify Api Tokens.

    get:
    Return details about a specific token.

    put:
    This can be used to activate a specific token by sending `is_activated: true`.
    Setting `is_activated` to true indicates agreement with refine.bio's
    [Terms of Use](https://www.refine.bio/terms) and
    [Privacy Policy](https://www.refine.bio/privacy).
    """

    model = APIToken
    lookup_field = "id"
    queryset = APIToken.objects.all()
    serializer_class = APITokenSerializer
