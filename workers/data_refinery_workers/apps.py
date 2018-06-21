from django.apps import AppConfig

class WorkersConfig(AppConfig):
    name = 'data_refinery_workers'
    verbose_name = "Data Refinery Workers"
    def ready(self):
        print("Ayyyy this is cool.")
