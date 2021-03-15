from rest_framework.pagination import LimitOffsetPagination


# Limit the max size of requests.
class LimitedLimitOffsetPagination(LimitOffsetPagination):
    max_limit = 1000
