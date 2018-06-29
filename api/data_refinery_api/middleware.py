import logging
from django.utils.deprecation import MiddlewareMixin

class SentryCatchBadRequestMiddleware(MiddlewareMixin):
    def process_response(self, request, response):
        if response.status_code < 400:
            return response

        from raven.contrib.django.models import client

        if not client.is_enabled():
            return response

        data = client.get_data_from_request(request)
        data.update({
            "level": logging.WARN,
            "logger": "BadRequestCatcher?",
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
