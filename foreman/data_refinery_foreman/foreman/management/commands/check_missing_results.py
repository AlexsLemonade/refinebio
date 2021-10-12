from django.core.management.base import BaseCommand

from data_refinery_common.models import Sample
from data_refinery_common.performant_pagination.pagination import PAGE_SIZE, PerformantPaginator


class Command(BaseCommand):
    def handle(self, *args, **options):
        samples = Sample.processed_objects.all()
        paginator = PerformantPaginator(samples, PAGE_SIZE)
        page = paginator.page()
        counter = 0
        while True:
            for sample in page.object_list:
                counter += 1
                if sample.results.count() == 0:
                    print(sample.accession_code)
            if not page.has_next():
                break
            else:
                page = paginator.page(page.next_page_number())

            if counter % 10000 == 0:
                print("Checked another 10000k samples.")
