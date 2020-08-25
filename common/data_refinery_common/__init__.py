from django.apps import AppConfig


class CommonAppConfig(AppConfig):
    name = 'data_refinery_common'

    # This will be called multiple times
    def ready(self):
        import data_refinery_common.signals

default_app_config = 'data_refinery_common.CommonAppConfig'