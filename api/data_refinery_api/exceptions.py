from rest_framework.exceptions import APIException
from rest_framework.serializers import ValidationError
from rest_framework.views import exception_handler


def custom_exception_handler(exc, context):
    data = {"error_type": "", "message": "", "details": ""}

    response = exception_handler(exc, context)

    if response is not None:
        if isinstance(exc, (InvalidFilters, InvalidData)):
            data["error_type"] = exc.default_code
            data["message"] = exc.default_detail
            data["details"] = exc.invalid_filters

    response.data = data
    return response


class InvalidFilters(APIException):
    status_code = 400
    default_detail = "You have supplied invalid filters. See `details` for a full list."
    default_code = "invalid_filters"
    invalid_filters = []

    def __init__(self, invalid_filters):
        super().__init__()
        self.invalid_filters = invalid_filters


class InvalidData(APIException):
    status_code = 400
    default_detail = ""
    default_code = "invalid_data"
    details = []

    def __init__(self, message=None, details=None):
        super().__init__(detail=message)
        self.details = details
