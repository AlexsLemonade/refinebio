import datetime

from django.conf import settings
from django.core.management.base import BaseCommand
from django.utils import timezone

import requests

from data_refinery_common.models import DatasetAnnotation


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
        days = options["days"]
        start_time = timezone.now() - datetime.timedelta(days=days)

        annotation_queryset = DatasetAnnotation.objects.filter(
            created_at__gt=start_time
        ).prefetch_related("dataset")
        annotations = [
            annotation
            for annotation in annotation_queryset
            if annotation.data["start"] and should_display_email(annotation.dataset.email_address)
        ]

        # Save the locations permanently, since IP addresses can cycle over time
        location_cache = dict()
        for annotation in annotation_queryset:
            if not "location" in annotation.data:
                ip = annotation.data["ip"]
                if not ip in location_cache:
                    location_cache[ip] = get_ip_location(ip)
                annotation.data["location"] = location_cache[ip]
                annotation.save()

        unique_users = list(set(annotation.dataset.email_address for annotation in annotations))
        unique_locations = list(set(annotation.data["location"] for annotation in annotations))

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
            is_returning_user = DatasetAnnotation.objects.filter(
                created_at__lt=start_time, dataset__email_address=email
            )
            text = "{0} | {1} downloads from {2}\n".format(email, downloads, locations)
            if is_returning_user:
                returning_users += text
            else:
                new_users += text

        if total_downloads > 0:
            fallback_text = "In the last {0} days, {1} users downloaded {2} datasets from {3} locations.".format(
                days, len(unique_users), total_downloads, len(unique_locations)
            )
        else:
            fallback_text = "There were no downloads in the last {0} days.".format(days)

        blocks = [
            {
                "type": "section",
                "text": {"type": "plain_text", "emoji": True, "text": fallback_text},
            }
        ]
        if new_users:
            blocks.append(
                {
                    "type": "section",
                    "text": {"type": "mrkdwn", "text": "*New users* \n" + new_users,},
                }
            )
        if returning_users:
            blocks.append(
                {
                    "type": "section",
                    "text": {"type": "mrkdwn", "text": "*Returning users* \n" + returning_users,},
                }
            )

        # Post to slack
        requests.post(
            settings.ENGAGEMENTBOT_WEBHOOK,
            json={
                "username": "EngagementBot",
                "icon_emoji": ":halal:",
                "channel": "#" + options["channel"],
                "text": fallback_text,
                "blocks": blocks,
            },
            headers={"Content-Type": "application/json"},
            timeout=10,
        )


def should_display_email(email: str) -> bool:
    """ Returns true if the given email is not associated with the CCDL suers """
    if not email:
        return False
    return not (
        email.startswith("cansav09")
        or email.startswith("arielsvn")
        or email.startswith("jaclyn.n.taroni")
        or email.startswith("kurt.wheeler")
        or email.startswith("greenescientist")
        or email.startswith("miserlou")
        or email.startswith("d.prasad")
        or email.endswith("@alexslemonade.org")
        or email is ("daniel.himmelstein@gmail.com")
        or email is ("dv.prasad991@gmail.com")
    )


def get_ip_location(remote_ip):
    try:
        data = requests.get("https://ipapi.co/" + remote_ip + "/json/", timeout=10).json()
        return "{0}, {1}".format(data["city"], data["country_name"])
    except Exception:
        return remote_ip
