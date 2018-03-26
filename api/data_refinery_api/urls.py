from django.conf.urls import url, include
from rest_framework.documentation import include_docs_urls
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.urlpatterns import format_suffix_patterns

from data_refinery_api.views import (
    ExperimentList, 
    ExperimentDetail,
    SampleList, 
    SampleDetail,
    OrganismList,
    PlatformList,
    InstitutionList
)


# This class provides a friendlier root API page.
class APIRoot(APIView):
    """
    Refine.bio API

    This open API provides access to all of the Refine.bio experiment, sample and result information.

    """
    def get(self, request):
        return Response({
            'experiments': reverse('experiments', request=request),
            'samples': reverse('samples', request=request),
            'organisms': reverse('organisms', request=request),
            'platforms': reverse('platforms', request=request),
            'institutions': reverse('institutions', request=request),
        })


urlpatterns = [
    # Endpoints / Self-documentation
    url(r'^experiments/$', ExperimentList.as_view(), name="experiments"),
    url(r'^experiments/(?P<pk>[0-9]+)/$', ExperimentDetail.as_view()),
    url(r'^samples/$', SampleList.as_view(), name="samples"),
    url(r'^samples/(?P<pk>[0-9]+)/$', SampleDetail.as_view()),
    url(r'^organisms/$', OrganismList.as_view(), name="organisms"),
    url(r'^platforms/$', PlatformList.as_view(), name="platforms"),
    url(r'^institutions/$', InstitutionList.as_view(), name="institutions"),

    url(r'^', APIRoot.as_view()),

    # Core API schema docs
    url(r'^docs/', include_docs_urls(title='Refine.bio API'))
]

urlpatterns = format_suffix_patterns(urlpatterns)
