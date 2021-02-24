import datetime

from django.conf import settings
from django.core.management.base import BaseCommand
from django.utils import timezone

import requests

from data_refinery_common.models import Dataset


def post_downloads_summary(days, channel):
    start_time = timezone.now() - datetime.timedelta(days=days)

    datasets = Dataset.processed_filtered_objects.filter(
        created_at__gt=start_time
    ).prefetch_related("datasetannotation_set")

    annotation_queryset = datasets.datasetannotation_set.all()

    # Save the locations permanently, since IP addresses can cycle over time
    location_cache = dict()
    for annotation in annotation_queryset:
        if "location" not in annotation.data:
            ip = annotation.data["ip"]
            if ip not in location_cache:
                location_cache[ip] = get_ip_location(ip)
            annotation.data["location"] = location_cache[ip]
            annotation.save()

    unique_users = list(set(dataset.email_address for dataset in datasets))
    unique_locations = list(set(annotation.data["location"] for annotation in annotation_queryset))

    new_users = ""
    returning_users = ""
    total_downloads = 0
    for email in unique_users:
        user_annotations = annotation_queryset.filter(dataset__email_address=email)
        downloads = user_annotations.count()
        total_downloads += downloads
        locations = ", ".join(
            list(set(annotation.data["location"] for annotation in user_annotations))
        )
        is_returning_user = Dataset.processed_filtered_objects.filter(
            created_at__lt=start_time, email_address=email
        )
        text = "{0} | {1} downloads from {2}\n".format(email, downloads, locations)
        if is_returning_user:
            returning_users += text
        else:
            new_users += text

    if total_downloads > 0:
        fallback_text = (
            f"In the last {days} days, {len(unique_users)} users downloaded {total_downloads}"
            f" datasets from {len(unique_locations)} locations."
        )
    else:
        fallback_text = f"There were no downloads in the last {days} days."

    blocks = [
        {"type": "section", "text": {"type": "plain_text", "emoji": True, "text": fallback_text}}
    ]
    if new_users:
        blocks.append(
            {"type": "section", "text": {"type": "mrkdwn", "text": "*New users* \n" + new_users}}
        )
    if returning_users:
        blocks.append(
            {
                "type": "section",
                "text": {"type": "mrkdwn", "text": "*Returning users* \n" + returning_users},
            }
        )

    # Post to slack
    requests.post(
        settings.ENGAGEMENTBOT_WEBHOOK,
        json={
            "username": "EngagementBot",
            "icon_emoji": ":halal:",
            "channel": "#" + channel,
            "text": fallback_text,
            "blocks": blocks,
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )


class Command(BaseCommand):
    help = "Post downloads summary to slack"

    def add_arguments(self, parser):
        parser.add_argument(
            "--days",
            type=int,
            default=7,  # default to a week
            help=("Number of days in the past for which to build the stats"),
        )
        parser.add_argument(
            "--channel",
            type=str,
            default="ccdl-general",
            help=("Optional parameter to choose the channel where the message will be posted."),
        )

    def handle(self, *args, **options):
        post_downloads_summary(options["days"], options["channel"])


def get_ip_location(remote_ip):
    try:
        data = requests.get("https://ipapi.co/" + remote_ip + "/json/", timeout=10).json()
        return "{0}, {1}".format(data["city"], data["country_name"])
    except (requests.exceptions.RequestException, ValueError, KeyError):
        return remote_ip
