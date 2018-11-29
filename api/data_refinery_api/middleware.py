import logging
from django.utils.deprecation import MiddlewareMixin
from data_refinery_common.utils import get_env_variable_gracefully

class SentryCatchBadRequestMiddleware(MiddlewareMixin):
    def process_response(self, request, response):
        if response.status_code < 400 or response.status_code == 404:
            return response

        from raven.contrib.django.models import client

        if not client.is_enabled():
            return response

        data = client.get_data_from_request(request)
        data.update({
            "level": logging.WARN,
            "logger": "BadRequestMiddleware",
        })
        message_template = "{status_code} code returned for URL: {url} {content}"

        content = ""
        try:
            content = "with message: {}".format(str(response.content))
        except:
            pass

        message = message_template.format(status_code=response.status_code,
                                          url=request.build_absolute_uri(),
                                          content=content)
        result = client.captureMessage(message=message, data=data)

        return response

    def process_exception(self, request, exception):
        from raven.contrib.django.models import client

        if not client.is_enabled():
            return

        data = client.get_data_from_request(request)
        data.update({
            "level": logging.WARN,
            "logger": "BadRequestMiddleware",
        })

        client.captureMessage(message=str(exception), data=data)

class RevisionMiddleware(MiddlewareMixin):
    def process_response(self, request, response):
        response['X-Source-Revision'] = get_env_variable_gracefully("SYSTEM_VERSION")
        return response
