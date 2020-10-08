def check_filters(view, special_filters=None):
    valid_filters = []

    if hasattr(view, "filterset_fields"):
        if view.filterset_fields:
            valid_filters += list(view.filterset_fields)

    if hasattr(view, "paginator"):
        if view.paginator:
            valid_filters += ["offset", "limit"]

    if hasattr(view, "ordering"):
        if view.ordering:
            valid_filters.append("ordering")

    if special_filters:
        valid_filters += special_filters

    invalid_filters = []
    for param in view.request.query_params:
        if param not in valid_filters:
            invalid_filters.append(param)

    return invalid_filters
