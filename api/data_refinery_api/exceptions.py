import pdb

from rest_framework.exceptions import APIException
from rest_framework.serializers import ValidationError
from rest_framework.views import exception_handler


def custom_exception_handler(exc, context):
    response = exception_handler(exc, context)

    if response is not None:

        data = {"error_type": "", "message": "", "details": ""}

        if isinstance(exc, (InvalidFilters, InvalidData, BadRequest)):
            data["error_type"] = exc.default_code
            data["message"] = exc.detail
            data["details"] = exc.details

            response.data = data

        if isinstance(exc, ValidationError):
            detail = exc.detail
            if isinstance(detail, dict):
                errors = []

                for key, value in detail.items():
                    value = value[0]
                    errors.append({"error_type": value.code, "message": value, "details": key})

                if len(errors) > 1:
                    data["error_type"] = "multiple_errors"
                    data[
                        "message"
                    ] = "Multiple errors have occurred. See `details` for a full list."
                    data["details"] = errors
                else:
                    data = errors[0]

                response.data = data

    return response


class InvalidFilters(APIException):
    status_code = 400
    default_detail = "You have supplied invalid filters. See `details` for a full list."
    default_code = "invalid_filters"
    details = []

    def __init__(self, message=None, invalid_filters=None):
        super().__init__(detail=message)
        self.details = invalid_filters


class InvalidData(APIException):
    status_code = 400
    default_detail = "The data that you have supplied is invalid"
    default_code = "invalid_data"
    details = []

    def __init__(self, message=None, details=None):
        super().__init__(detail=message)
        self.details = details


class BadRequest(APIException):
    status_code = 400
    default_detail = ""
    default_code = "bad_request"
    details = ""

    def __init__(self, message=None, error_type=None):
        super().__init__(detail=message)
        if error_type:
            self.default_code = error_type
