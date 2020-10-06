from rest_framework.exceptions import APIException


class InvalidFilters(APIException):
    status_code = 400
    default_detail = "You have supplied invalid filters"
    default_code = "invalid_query_parameters"
