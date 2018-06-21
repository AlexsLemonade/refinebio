from django.apps import AppConfig

class ApiConfig(AppConfig):
    name = 'data_refinery_api'
    verbose_name = "Data Refinery API"
    def ready(self):
        print("Ayyyy this is cool.")
