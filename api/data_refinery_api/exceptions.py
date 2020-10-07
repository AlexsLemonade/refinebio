from rest_framework.exceptions import APIException
from rest_framework.views import exception_handler


def custom_exception_handler(exc, context):
    response = exception_handler(exc, context)

    if response:
        if isinstance(exc, APIException):
            response.data["error_code"] = exc.detail.code

        if isinstance(exc, InvalidFilters):
            response.data["invalid_filters"] = exc.invalid_filters

    return response


class InvalidFilters(APIException):
    status_code = 400
    default_detail = "You have supplied invalid filters. A list of the filters that were invalid are provided in the invalid_filters key."
    default_code = "invalid_filters"

    def __init__(self, invalid_filters):
        super().__init__()
        self.invalid_filters = invalid_filters
