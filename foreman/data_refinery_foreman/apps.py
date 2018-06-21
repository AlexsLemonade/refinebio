from django.apps import AppConfig
from data_refinery_common.utils import get_env_variable
from raven import Client
import logging

class ForemanConfig(AppConfig):
    name = 'data_refinery_foreman'
    verbose_name = "Data Refinery Foreman"
    def ready(self):
        raven_dsn = get_env_variable("RAVEN_DSN", "")
        if raven_dsn == "":
            # If there is no DSN for Raven, it will complain every
            # time it should have sent an event.
            raven_logger = logging.getLogger('raven.base.Client')
            raven_logger.setLevel(logging.CRITICAL)
