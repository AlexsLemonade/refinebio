import datetime

from django.db.models import Count
from django.db.models.functions import ExtractWeek
from django.utils import timezone

from data_refinery_common.models import Dataset

START_DATE = datetime.datetime(2021, 1, 1, tzinfo=timezone.utc)

qs = (
    Dataset.processed_filtered_objects.filter(created_at__gt=START_DATE)
    .annotate(
        week=ExtractWeek("created_at"),
    )
    .values("week")
    .annotate(n=Count("pk"))
    .order_by("week")
)

print("week_number, downloads")
for week in qs:
    print(f"{week['week']}, {week['n']}")

total_for_year = Dataset.processed_filtered_objects.filter(created_at__gt=START_DATE).count()

print("Total for year: {total_for_year}")
