from django.apps import AppConfig

class ApiConfig(AppConfig):
    name = 'data_refinery_api'
    verbose_name = "Data Refinery API"
    def ready(self):
        raven_dsn = get_env_variable("RAVEN_DSN", "")
        if raven_dsn == "":
            # If there is no DSN for Raven, it will complain every
            # time it should have sent an event.
            raven_logger = logging.getLogger('raven.base.Client')
            raven_logger.setLevel(logging.CRITICAL)
