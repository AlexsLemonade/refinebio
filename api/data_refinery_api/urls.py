from django.conf.urls import url
from rest_framework.documentation import include_docs_urls
from rest_framework.urlpatterns import format_suffix_patterns
from data_refinery_api.views import ExperimentList, ExperimentDetail

urlpatterns = [
	# Endpoints / Self-documentation
    url(r'^experiments/$', ExperimentList.as_view()),
    url(r'^experiments/(?P<pk>[0-9]+)/$', ExperimentDetail.as_view()),

    # Core API schema docs
	url(r'^docs/', include_docs_urls(title='Refine.bio API'))
]

urlpatterns = format_suffix_patterns(urlpatterns)